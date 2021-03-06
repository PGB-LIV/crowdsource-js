/**
 * Copyright 2019 University of Liverpool
 * 
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not
 * use this file except in compliance with the License. You may obtain a copy of
 * the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
 * License for the specific language governing permissions and limitations under
 * the License.
 */

/**
 * Web Worker
 */
this.search = new MsSearch();
this.onmessage = function(event) {
	"use strict";
	var result = this.search.search(event.data);

	if (typeof WorkerGlobalScope !== 'undefined'
			&& self instanceof WorkerGlobalScope) {
		postMessage(result);
		return;
	}

	sendResult(result);
};

function MsSearch(data) {
	"use strict";

	var VERSION = 'INSERT_BUILD_VERSION';
	var H2O_MASS = 1.00782503223 + 1.00782503223 + 15.99491461957;
	var CO_MASS = 12 + 15.99491461957;
	var NH3_MASS = 14.00307400443 + 1.00782503223 + 1.00782503223 + 1.00782503223;
	var PHOS_LOSS_MASS = 97.97689557339;
	var PROTON_MASS = 1.007276466879;
	var AA_MASS = {
		A : 71.0371137852,
		R : 156.1011110241,
		N : 114.0429274414,
		D : 115.0269430243,
		C : 103.0091849596,
		E : 129.0425930888,
		Q : 128.0585775058,
		G : 57.0214637207,
		H : 137.0589118585,
		I : 113.0840639785,
		L : 113.0840639785,
		K : 128.0949630152,
		M : 131.0404850885,
		F : 147.0684139141,
		P : 97.0527638496,
		S : 87.0320284047,
		T : 101.0476784692,
		W : 186.0793129507,
		Y : 163.0633285336,
		V : 99.0684139141,
		U : 150.9536355852
	};

	this.workUnit = null;
	this.spectra = null;
	this.probability = null;
	this.getTolerance = null;

	this.initialise = function(data) {
		this.workUnit = data;
		this.spectra = this.indexSpectra(this.workUnit.fragments);
		this.probability = this.getProbability(this.workUnit.fragments);

		if (this.workUnit.fragTolUnit === 'ppm') {
			this.getTolerance = function(mass) {
				return mass * 0.000001 * this.workUnit.fragTol;
			};
		} else {
			this.getTolerance = function(mass) {
				return this.workUnit.fragTol;
			};
		}
	};

	this.search = function(data) {
		this.initialise(data);

		var searchStart = Date.now();
		var peptideScore = [];

		var searchResult;
		for (var peptideIndex = 0; peptideIndex < this.workUnit.peptides.length; peptideIndex++) {
			searchResult = this
					.searchSequence(this.workUnit.peptides[peptideIndex]);
			peptideScore[peptideIndex] = searchResult;
		}

		peptideScore.sort(function(a, b) {
			return b.score - a.score;
		});

		var resultObject = {
			uid : this.workUnit.uid,
			processTime : Date.now() - searchStart,
			peptides : []
		};

		for (var r = 0; r < this.workUnit.peptides.length; r++) {
			resultObject.peptides[r] = peptideScore[r];
		}

		return resultObject;
	};

	this.searchSequence = function(peptide) {
		var totalModNum = this.getTotalModNum(peptide.mods);

		// Array of ModLoc objects (all possible locations of all modifications
		// reported) ModLoc objects are {possLoc,modIndex,vModMass}]
		var modLocs = this.getAllModLocs(peptide);

		// flushes out resObj arrays to correct size
		var resObj = this.getInitialResObj(peptide);

		// just check..
		if (modLocs.length < totalModNum) {
			totalModNum = modLocs.length;
		}

		var subIonsets = this.getSequenceIonSets(peptide, modLocs, totalModNum);
		var subIonMatches = [];

		for (var s = 0; s < subIonsets.length; s++) {
			// for each ionset we log matches with ms2 fragments
			subIonsets[s] = this.matchSpectraWithIonSets(subIonsets[s],
					peptide.sequence, this.workUnit.z);

			// score the ionset (matched ions and ion ladders)
			subIonMatches[s] = this.getMatchCount(subIonsets[s]);
		}

		var bestScore = this.getBestScore(peptide.sequence, subIonsets,
				subIonMatches);

		// currScoreObj = best score, ion-match and modPos of this best
		// score
		resObj.S = Math.round(bestScore.score * 100000) / 100000;
		resObj.IM = bestScore.ionsMatchSum;

		// sort out resobj mod positions
		for (var i = 0; i < bestScore.modPos.length; i++) {
			// currScoreObj.modPos is an array of ModLocs
			var m = bestScore.modPos[i].modIndex;
			var p = bestScore.modPos[i].possLoc;
			resObj.mods[m].P.push(p); // to resObj.mod[].position;
		}

		// For the moment I have to bulk out the position field to be equal
		// to number of mods (if >=200 then I haven't found/cannot confirm a
		// position)
		for (var index = 0; index < resObj.mods.length; index++) {
			if (resObj.mods[index].P.length < peptide.mods[index].num) {
				var n = peptide.mods[index].num - resObj.mods[index].P.length;
				for (var t = 0; t < n; t++) {
					resObj.mods[index].P.push(200 + t);
				}
			}
		}

		return resObj;
	};

	this.getBestScore = function(sequence, subIonsets, subIonMatches) {
		var bestScore = {
			modPos : [],
			ionsMatched : 0,
			ionsMatchSum : 0,
			score : 0
		};

		var currScore;
		for (var setIndex = 0; setIndex < subIonsets.length; setIndex++) {
			var aCount = subIonsets[setIndex].aIons.length - 1;
			var ionCount = ((sequence.length - 2) * 2) + aCount;

			currScore = 0;
			if (subIonMatches[setIndex] > 0) {
				currScore = this.BinomialP(1 - this.probability, ionCount,
						ionCount - subIonMatches[setIndex] - 1);

				currScore = currScore === Infinity ? 0 : -10
						* Math.log10(currScore);
			}

			// Select best candidate
			if (currScore < bestScore.score
					|| (currScore === bestScore.score && this
							.getMatchSum(subIonsets[setIndex]) < bestScore.ionsMatchSum)) {
				continue;
			}

			bestScore.score = currScore;
			bestScore.modPos = subIonsets[setIndex].modPos.slice();
			// place modPos structure [ModLoc] is kept with score
			bestScore.ionsMatched = subIonMatches[setIndex];
			bestScore.ionsMatchSum = this.getMatchSum(subIonsets[setIndex]);
		}

		return bestScore;
	};

	// common
	this.matchSpectraWithIonSets = function(ionSet, sequence, maxCharge) {
		this.matchSpectraWithIonSet(ionSet.aIons, sequence, maxCharge);
		this.matchSpectraWithIonSet(ionSet.bIons, sequence, maxCharge);
		this.matchSpectraWithIonSet(ionSet.yIons, sequence, maxCharge);

		return ionSet;
	};

	this.matchSpectraWithIonSet = function(ions, sequence, maxCharge) {
		var possibleIons;
		var losses = [];
		var neutralIon;

		for (var i = 0; i < ions.length; i++) {
			neutralIon = ions[i];

			this.updateLosses(losses, neutralIon, sequence.charAt(i));

			if (i === 0) {
				continue;
			}

			possibleIons = this.getIons(neutralIon, maxCharge, losses);

			this.matchSpectraWithPossibleIons(possibleIons, neutralIon);
		}
	};

	this.matchSpectraWithPossibleIons = function(possibleIons, neutralIon) {
		var possibleIon;
		for (var ionIndex = 0; ionIndex < possibleIons.length; ionIndex++) {
			possibleIon = possibleIons[ionIndex];

			if (!this.isMassInSpectra(possibleIon.mass)) {
				continue;
			}

			if (possibleIon.baseMatch === 1) {
				neutralIon.baseMatch = 1;
			}

			neutralIon.match++;
		}
	};

	this.updateLosses = function(losses, neutralIon, residue) {
		if (neutralIon.modification === 21) {
			losses.Phos = true;
		}

		if (losses.H2O !== true
				&& (residue === 'S' || residue === 'T' || residue === 'E' || residue === 'D')) {
			losses.H2O = true;
		}

		if (losses.NH3 !== true
				&& (residue === 'R' || residue === 'K' || residue === 'N' || residue === 'Q')) {
			losses.NH3 = true;
		}
	};

	/**
	 * Is ion mass (mass) found in ms2 spectrum (spectra) returns zero if no
	 * match else first match intensity. also store deltaM (how close to
	 * precision window we are) common
	 */
	this.isMassInSpectra = function(mass) {
		var tolerance = this.getTolerance(mass);

		var idMax = Math.floor(mass + tolerance);
		for (var id = Math.floor(mass - tolerance); id <= idMax; id++) {
			if (typeof this.spectra[id] === 'undefined'
					|| !this.isInSpectraWindow(this.spectra[id], mass,
							tolerance)) {
				continue;
			}

			return true;
		}

		return false;
	};

	this.isInSpectraWindow = function(spectra, mass, precision) {
		var deltaM = 0;
		for (var m = 0; m < spectra.length; m++) {
			deltaM = Math.abs(mass - spectra[m].mz);

			if (deltaM < precision) {
				return true;
			}
		}

		return false;
	};

	this.checkforFixedPTM = function(residue) {
		var massShift = 0;
		for (var m = 0; m < this.workUnit.fixedMods.length; m++) {
			if (this.workUnit.fixedMods[m].residues.length > 1) {
				// just to check...
				console.log('Fixed Mod needs looking at '
						+ this.workUnit.fixedMods[m].residues);
			}

			if (this.workUnit.fixedMods[m].residues === residue) {
				massShift += this.workUnit.fixedMods[m].mass;
			}
		}

		return massShift;
	};

	this.getMatchCountIons = function(ions) {
		var count = 0;
		for (var i = 0; i < ions.length; i++) {
			if (ions[i].baseMatch > 0) {
				count++;
			} else if (ions[i].match > 0) {
				count += 0.45;
			}
		}

		return count;
	};

	this.getMatchCount = function(ionSet) {
		// A ions not counted
		var bCount = this.getMatchCountIons(ionSet.bIons);
		var yCount = this.getMatchCountIons(ionSet.yIons);

		return bCount + yCount;
	};

	this.getMatchSumIons = function(ions) {
		var count = 0;
		for (var i = 0; i < ions.length; i++) {
			count += ions[i].match;
		}

		return count;
	};

	this.getMatchSum = function(ionSet) {
		var aCount = this.getMatchSumIons(ionSet.aIons);
		var bCount = this.getMatchSumIons(ionSet.bIons);
		var yCount = this.getMatchSumIons(ionSet.yIons);

		return aCount + bCount + yCount;
	};

	/**
	 * Gets the total modifiable sites from the peptide mod list
	 */
	this.getTotalModNum = function(mods) {
		var mnum = 0;
		for (var mod = 0; mod < mods.length; mod++) {
			mnum += mods[mod].num;
		}

		return mnum;
	};

	this.getInitialResObj = function(peptide) {
		var resObj = {
			id : peptide.id,
			IM : 0,
			S : 0,
			mods : []
		};

		if (typeof peptide.mods === 'undefined') {
			return resObj;
		}

		for (var mod = 0; mod < peptide.mods.length; mod++) {
			var modp = {
				id : peptide.mods[mod].id,
				P : []
			};

			resObj.mods[mod] = modp;
		}

		return resObj;
	};

	this.getAllModLocs = function(peptide) {
		var mlocs = [];
		if (typeof peptide.mods === 'undefined') {
			return mlocs;
		}

		for (var mod = 0; mod < peptide.mods.length; mod++) {
			var possLocs = this.createIndexofPossibleLocations3(
					peptide.sequence, peptide.mods[mod].residues);

			for (var i = 0; i < possLocs.length; i++) {
				var modObj = new ModLoc(possLocs[i], mod,
						peptide.mods[mod].mass, peptide.mods[mod].id);
				mlocs.push(modObj);
			}
		}

		mlocs.sort(function(a, b) {
			return a.possLoc - b.possLoc;
		});

		return mlocs;
	};

	this.getSequenceIonSets = function(modifiedSequence, modlocs, num) {
		var ionSetArray = [];

		// Unmodified
		if (num === 0) {
			ionSetArray[0] = this.getSequenceIons(modifiedSequence, []);

			return ionSetArray;
		}

		// Single modification
		if (num === 1) {
			for (var ml = 0; ml < modlocs.length; ml++) {
				ionSetArray[ml] = this.getSequenceIons(modifiedSequence,
						[ modlocs[ml] ]);
			}

			return ionSetArray;
		}

		if (modlocs.length < num) {
			return ionSetArray;
		}

		var combArray = this.createArrayPossibleCombinations3(modlocs,
				modlocs.length, num); // combarray
		// now in form of [modlocs]

		for (var c = 0; c < combArray.length; c++) {
			var mlocs = combArray[c]; // mlocs = confirmed mods + combination
			// to

			var modFreq = [];
			for (var i = 0; i < modifiedSequence.mods.length; i++) {
				modFreq[i] = 0;
			}

			for (i = 0; i < mlocs.length; i++) {
				modFreq[mlocs[i].modIndex]++;
			}

			var isMatch = true;
			for (i = 0; i < modifiedSequence.mods.length; i++) {
				if (modFreq[i] !== modifiedSequence.mods[i].num) {
					isMatch = false;
				}
			}

			if (!isMatch) {
				continue;
			}

			var ions = this.getSequenceIons(modifiedSequence, mlocs);

			ionSetArray.push(ions);
		}

		return ionSetArray;
	};

	this.getChargedMass = function(mass, charge) {
		// Add proton
		var chargedMass = mass + (PROTON_MASS * charge);

		return chargedMass / charge;
	};

	this.getIonsA = function(sequence, modificationLocations) {
		var neutralMass = this.checkforFixedPTM('[');
		neutralMass += CO_MASS;

		var residue;
		var ions = [];
		var ionObj;
		for (var a = 0; a < sequence.length; a++) {
			ionObj = new Ion();
			residue = sequence.charAt(a);
			neutralMass += AA_MASS[residue] + this.checkforFixedPTM(residue);

			if (neutralMass > 800) {
				break;
			}

			for (var loopIndex = 0; loopIndex < modificationLocations.length; loopIndex++) {
				if (a === (modificationLocations[loopIndex].possLoc - 1)) {
					neutralMass += modificationLocations[loopIndex].vModMass;
					ionObj.modification = modificationLocations[loopIndex].id;
				}
			}

			ionObj.mass = neutralMass;
			ions[a] = ionObj;
		}

		return ions;
	};

	this.getIonsB = function(sequence, modificationLocations) {
		var neutralMass = this.checkforFixedPTM('[');

		var residue;
		var ions = [];
		var ionObj;
		for (var b = 0; b < sequence.length - 1; b++) {
			ionObj = new Ion();
			residue = sequence.charAt(b);
			neutralMass += AA_MASS[residue] + this.checkforFixedPTM(residue);

			for (var loopIndex = 0; loopIndex < modificationLocations.length; loopIndex++) {
				if (b === (modificationLocations[loopIndex].possLoc - 1)) {
					neutralMass += modificationLocations[loopIndex].vModMass;
					ionObj.modification = modificationLocations[loopIndex].id;
				}
			}

			ionObj.mass = neutralMass;
			ions[b] = ionObj;
		}

		return ions;
	};

	this.getIonsY = function(sequence, modificationLocations) {
		var neutralMass = H2O_MASS + this.checkforFixedPTM(']');

		var residue;
		var ions = [];
		var ionObj;
		for (var y = sequence.length - 1; y > 0; y--) {
			ionObj = new Ion();
			residue = sequence.charAt(y);
			neutralMass += AA_MASS[residue] + this.checkforFixedPTM(residue);

			for (var loopIndex = 0; loopIndex < modificationLocations.length; loopIndex++) {
				if (y === (modificationLocations[loopIndex].possLoc - 1)) {
					neutralMass += modificationLocations[loopIndex].vModMass;
					ionObj.modification = modificationLocations[loopIndex].id;
				}
			}

			ionObj.mass = neutralMass;
			ions.push(ionObj);
		}

		return ions;
	};

	this.getIons = function(neutralIon, maxCharge, losses) {
		var ions = [];
		var ionObj;
		for (var charge = 1; charge <= maxCharge; charge++) {
			ionObj = new Ion();
			ionObj.mass = this.getChargedMass(neutralIon.mass, charge);
			ionObj.baseMatch = 1;
			ions.push(ionObj);

			// Water
			if (losses.H2O === true) {
				ionObj = new Ion();
				ionObj.mass = this.getChargedMass(neutralIon.mass - H2O_MASS,
						charge);
				ions.push(ionObj);
			}

			// Amonia
			if (losses.NH3 === true) {
				ionObj = new Ion();
				ionObj.mass = this.getChargedMass(neutralIon.mass - NH3_MASS,
						charge);
				ions.push(ionObj);
			}

			// Water + Amonia
			if (losses.H2O === true && losses.NH3 === true) {
				ionObj = new Ion();
				ionObj.mass = this.getChargedMass(neutralIon.mass
						- (H2O_MASS + NH3_MASS), charge);
				ions.push(ionObj);
			}

			if (losses.Phos === true) {
				ionObj = new Ion();
				ionObj.mass = this.getChargedMass(neutralIon.mass
						- PHOS_LOSS_MASS, charge);
				ions.push(ionObj);
			}
		}

		return ions;
	};

	this.getSequenceIons = function(modifiedSequence, mlocs) {
		var ionSet = new Ionset();

		ionSet.modPos = mlocs.slice();
		ionSet.aIons = this.getIonsA(modifiedSequence.sequence, mlocs);
		ionSet.bIons = this.getIonsB(modifiedSequence.sequence, mlocs);
		ionSet.yIons = this.getIonsY(modifiedSequence.sequence, mlocs);

		return ionSet;
	};

	this.indexSpectra = function(fragments) {
		var spectra = [];
		var key;
		for (var i = 0; i < fragments.length; i++) {
			key = Math.floor(fragments[i].mz);
			if (typeof spectra[key] === 'undefined') {
				spectra[key] = [];
			}

			spectra[key].push(fragments[i]);
		}

		return spectra;
	};

	// returns integer array of all possible locations of mod +1
	// NB position = sequence position + 1 ... 0 = '[' 1 to len = = to len-1 and
	// length+1 = ] ]
	// input sequence (eg 'PEPTIDE') and residues (eg 'CKV');
	this.createIndexofPossibleLocations3 = function(sequence, residues) {
		var posloc = []; // array of possible locations for this mod
		var pos = 0;
		var res = '';
		// special case for residues '[' annd ']'
		for (var r = 0; r < residues.length; r++) {
			res = residues[r];
			if (res === '[') {
				pos = 0;
				posloc.push(pos);
			} else if (res === ']') {
				pos = sequence.length + 1;
				posloc.push(pos);
			} else {
				for (pos = 0; pos < sequence.length; pos++) {
					if (sequence[pos] === res) {
						posloc.push(pos + 1); // using modified location 1 =
						// position zero
					}
				}
			}
		}

		return posloc;
	};

	// -----------------------------------------------------------------------------
	// Combination & Factorial code

	this.createArrayPossibleCombinations3 = function(locArray, n, k) {
		var result = [];
		var values = [];
		for (var i = 0; i <= n; i++) {
			values[i] = i;
		}

		// initialize permutations
		var perm = [];
		for (i = 0; i < n; i++) {
			if (i < k) {
				perm[i] = 1;
			} else {
				perm[i] = 0;
			}
		}

		perm.sort();

		whileloop: while (true) {
			// save subresult
			var subresult = [];
			for (i = 0; i < n; i++) {
				if (perm[i] === 1) {
					subresult.push(values[i]);
				}
			}

			result.push(subresult);
			// get next permutation
			for (i = n - 1; i > 0; i--) {
				if (perm[i - 1] === 1) {
					continue;
				}

				if (perm[i] === 1) {
					perm[i - 1] = 1;
					perm[i] = 0;
					perm = perm.slice(0, i).concat(perm.slice(i).sort());
					continue whileloop;
				}
			}
			// no additional permutations exist
			break whileloop;
		}

		// result = array of integers that correspond to index in location array
		var retObj = [];
		for (i = 0; i < result.length; i++) {
			retObj[i] = [];
			for (var j = 0; j < k; j++) {
				retObj[i][j] = locArray[result[i][j]]; // create actual mod
				// location.
			}
		}

		return retObj;
	};

	this.getCombinations = function(n, k) {
		var nfac = getFactorial(n);
		var kfac = getFactorial(k);
		var nsubkfac = getFactorial(n - k);

		return nfac / (kfac * nsubkfac);
	};

	this.getFactorial = function(f) {
		var r = 1;
		for (; f > 0; r *= f, f--) {
			// Logic in iterator
		}

		return r;
	};

	// SRC: http://www.ciphersbyritter.com/JAVASCRP/BINOMPOI.HTM
	this.BinomialP = function(p, n, k) {
		// term-by-term
		if ((k > n) || (p >= 1)) {
			return 1;
		}

		var q = 1 - p;
		var n1p = (n + 1) * p;

		var t = n * Math.log(q); // k = 0
		var r = Math.exp(t);
		var j = 1;
		while (j <= k) {
			t += Math.log(1 + (n1p - j) / (j * q));
			r += Math.exp(t);
			j++;
		}

		return r;
	};

	this.getProbability = function(fragments) {
		var N = fragments.length;

		var D;

		if (this.workUnit.fragTolUnit === 'ppm') {
			D = 0.000001 * this.workUnit.fragTol;
		} else {
			D = this.workUnit.fragTol;
		}

		var W = 0;
		for (var fragId = 0; fragId < N; fragId++) {
			W = Math.max(W, fragments[fragId].mz);
		}

		W = Math.ceil(W / 100) * 100;

		return N * D / W;
	};
}

function Ionset() {
	/**
	 * an array of ModLocs
	 */
	this.modPos = [];

	/**
	 * an array of Ion
	 */
	this.aIons = [];

	/**
	 * an array of Ion
	 */
	this.bIons = [];

	/**
	 * an array of Ion
	 */
	this.yIons = [];
}

function Ion() {
	this.mass = 0;
	this.baseMatch = 0;
	this.match = 0;
	this.intensity = 0;
	this.deltaM = 0;
	this.modification = 0;
}

/**
 * used in second/third phase
 */
function ModLoc(ploc, modindex, modmass, id) {
	/**
	 * location on peptide
	 */
	this.possLoc = ploc;

	/**
	 * modification index (within Work unit Peptide description)
	 */
	this.modIndex = modindex;

	/**
	 * mass of this modification (convenience)
	 */
	this.vModMass = modmass;

	/**
	 * ID of modification
	 */
	this.id = id;
}

if (typeof Math.log10 === 'undefined') {
	Object.prototype.log10 = function(n) {
		return Math.log(n) / Math.log(10);
	};
}
