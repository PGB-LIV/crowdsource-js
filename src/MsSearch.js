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
this.onmessage = function(event) {
	"use strict";
	var search = new MsSearch(event.data);
	var resultObject = search.search();

	if (typeof WorkerGlobalScope !== 'undefined'
			&& self instanceof WorkerGlobalScope) {
		postMessage(resultObject);
		return;
	}

	sendResult(resultObject);
};

function MsSearch(data) {
	"use strict";

	var H2O_MASS = 1.00782503223 + 1.00782503223 + 15.99491461957;
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

	this.workUnit = data;

	if (this.workUnit.fragTolUnit === 'ppm') {
		this.getTolerance = function(mass) {
			return mass * 0.000001 * this.workUnit.fragTol;
		};
	} else {
		this.getTolerance = function(mass) {
			return this.workUnit.fragTol;
		};
	}

	this.search = function() {
		var searchStart = Date.now();
		var allPeptideScores = [];
		var spectra = this.indexSpectra(this.workUnit.fragments);

		for (var peptideIndex = 0; peptideIndex < this.workUnit.peptides.length; peptideIndex++) {
			var result = this.searchSequence(spectra,
					this.workUnit.peptides[peptideIndex]);
			allPeptideScores[peptideIndex] = result;
		}

		allPeptideScores.sort(function(a, b) {
			return b.score - a.score;
		});

		var resultObject = {
			uid : this.workUnit.uid,
			processTime : Date.now() - searchStart,
			peptides : []
		};

		for (var r = 0; r < this.workUnit.peptides.length; r++) {
			resultObject.peptides[r] = allPeptideScores[r];
		}

		return resultObject;
	};

	this.searchSequence = function(spectra, currPeptide) {
		var currScoreObj = {
			modPos : [],
			ionsMatched : 0,
			ionsMatchSum : 0,
			score : 0
		};

		var totalModNum = this.getTotalModNum(currPeptide);

		// Array of ModLoc objects (all possible locations of all modifications
		// reported) ModLoc objects are {possLoc,modIndex,vModMass}]
		var modLocs = this.getAllModLocs(currPeptide);

		// flushes out resObj arrays to correct size
		var resObj = this.getInitialResObj(currPeptide);

		// just check..
		if (modLocs.length < totalModNum) {
			totalModNum = modLocs.length;
		}

		var subIonsets = this.getSequenceIonSets(currPeptide, modLocs,
				totalModNum);

		// an array of ionsets to score containing num number of mods
		var subIonMatches = [];

		var ptmRsP = this.getFactorP(this.workUnit.fragments);
		for (var s = 0; s < subIonsets.length; s++) {
			// for each ionset we log matches with ms2 fragments
			subIonsets[s] = this.matchSpectraWithIonSets(spectra,
					subIonsets[s], currPeptide.sequence, this.workUnit.z);

			// score the ionset (matched ions and ion ladders)
			subIonMatches[s] = this.getMatchCount(subIonsets[s]);
		}

		var ionCount = (currPeptide.sequence.length - 1) * 2;
		for (s = 0; s < subIonsets.length; s++) {
			var score = 0;
			if (subIonMatches[s] > 0) {
				score = this.BinomialP(1 - ptmRsP, ionCount, ionCount
						- subIonMatches[s] - 1);

				score = score === Infinity ? 0 : -10 * Math.log10(score);
			}

			// Select best candidate
			if (score > currScoreObj.score
					|| (score === currScoreObj.score && this
							.getMatchSum(subIonsets[s]) >= currScoreObj.ionsMatchSum)) {
				currScoreObj.score = score;
				currScoreObj.modPos = subIonsets[s].modPos.slice();
				// place modPos structure [ModLoc] is kept with score
				currScoreObj.ionsMatched = subIonMatches[s];
				currScoreObj.ionsMatchSum = this.getMatchSum(subIonsets[s]);
			}
		}

		// currScoreObj = best score, ion-match and modPos of this best
		// score
		resObj.S = Math.round(currScoreObj.score * 100000) / 100000;
		resObj.IM = currScoreObj.ionsMatchSum;

		// sort out resobj mod positions
		for (var i = 0; i < currScoreObj.modPos.length; i++) {
			// currScoreObj.modPos is an array of ModLocs
			var m = currScoreObj.modPos[i].modIndex;
			var p = currScoreObj.modPos[i].possLoc;
			resObj.mods[m].P.push(p); // to resObj.mod[].position;
		}

		// For the moment I have to bulk out the position field to be equal
		// to number of mods (if >=200 then I haven't found/cannot confirm a
		// position)
		for (var index = 0; index < resObj.mods.length; index++) {
			if (resObj.mods[index].P.length < currPeptide.mods[index].num) {
				var n = currPeptide.mods[index].num
						- resObj.mods[index].P.length;
				for (var t = 0; t < n; t++) {
					resObj.mods[index].P.push(200 + t);
				}
			}
		}

		return resObj;
	};

	// common
	this.matchSpectraWithIonSets = function(spectra, ionSet, sequence,
			maxCharge) {
		this.matchSpectraWithIonSet(spectra, ionSet.bIons, sequence, maxCharge);
		this.matchSpectraWithIonSet(spectra, ionSet.yIons, sequence, maxCharge);

		return ionSet;
	};

	this.matchSpectraWithIonSet = function(spectra, ions, sequence, maxCharge) {
		var possibleIons;
		var losses = [];
		var neutralIon;

		for (var i = 0; i < ions.length; i++) {
			neutralIon = ions[i];

			this.updateLosses(losses, neutralIon, sequence.charAt(i));
			possibleIons = this.getIons(neutralIon, maxCharge, losses);

			this
					.matchSpectraWithPossibleIons(spectra, possibleIons,
							neutralIon);
		}
	};

	this.matchSpectraWithPossibleIons = function(spectra, possibleIons,
			neutralIon) {
		var possibleIon;
		for (var ionIndex = 0; ionIndex < possibleIons.length; ionIndex++) {
			possibleIon = possibleIons[ionIndex];

			if (!this.isMassInSpectra(spectra, possibleIon.mass)) {
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
	this.isMassInSpectra = function(spectra, mass) {
		var tolerance = this.getTolerance(mass);

		var idMax = Math.floor(mass + tolerance);
		for (var id = Math.floor(mass - tolerance); id <= idMax; id++) {

			if (typeof spectra[id] === 'undefined'
					|| !this.isInSpectraWindow(spectra[id], mass, tolerance)) {
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

	this.getMatchCount = function(ionSet) {
		var bCount = 0;
		for (var b = 0; b < ionSet.bIons.length; b++) {
			if (ionSet.bIons[b].baseMatch > 0) {
				bCount++;
			} else if (ionSet.bIons[b].match > 0) {
				bCount += 0.75;
			}
		}

		var yCount = 0;
		for (var y = 0; y < ionSet.yIons.length; y++) {
			if (ionSet.yIons[y].baseMatch > 0) {
				yCount++;
			} else if (ionSet.yIons[y].match > 0) {
				yCount += 0.75;
			}
		}

		return bCount + yCount;
	};

	this.getMatchSum = function(ionSet) {
		var bCount = 0;
		for (var b = 0; b < ionSet.bIons.length; b++) {
			bCount += ionSet.bIons[b].match;
		}

		var yCount = 0;
		for (var y = 0; y < ionSet.yIons.length; y++) {
			yCount += ionSet.yIons[y].match;
		}

		return bCount + yCount;
	};

	/**
	 * Gets the total modifiable sites from the modification list
	 * 
	 * @param peptide
	 *            Object to read mods from
	 * @returns int
	 */
	this.getTotalModNum = function(peptide) {
		var mnum = 0;
		for (var mod = 0; mod < peptide.mods.length; mod++) {
			mnum += peptide.mods[mod].num;
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

		// why do we sort it? to shuffle the mods?
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
		// do a calculation first to see if do-able combinations rapidly get out
		// of hand.
		// number of combinations = n!/(k!(n-k)!)
		// if (getCombinations(modlocs.length, num) > MUCH_TOO_MUCH) {
		// console.log('Bailing: Array too large for me.');
		// return ionSetArray;
		// }

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

			// try
			// console.log('mlocs:'+JSON.stringify(mlocs)+'
			// len:'+modlocs.length+'
			// num:'+num+'conpos:'+conpos.length);
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

	this.getFactorP = function(fragments) {
		var ptmRsN = fragments.length;
		var ptmRsD = 0.5;
		var ptmRsW = 0;
		for (var fragId = 0; fragId < ptmRsN; fragId++) {
			ptmRsW = Math.max(ptmRsW, fragments[fragId].mz);
		}

		// TODO: Remove hard coded window size
		ptmRsW = Math.ceil(ptmRsW / 100) * 100;

		return ptmRsN * ptmRsD / ptmRsW;
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