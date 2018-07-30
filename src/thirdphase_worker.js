"use strict";
/**
 * Copyright 2016 University of Liverpool
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
var PHOSPHORYLATION = 21; // global equate. (may change dependent on
// Uniprot!!!!)
var PHOSPHO_MASS = 79.966331;
var PHOSPHO_LOSS = 97.966; // Phosphorylation is a special case hence special
// globals ....
var PHOSPHO2_LOSS = 115.987;
// in actual fact probably all modifications require special case scoring due to
// different proton affinities etc.
// but phophorylation definitely does so rather than have a 500 entry jump table
// /switch statement we'll try this first.

/**
 * 
 * @param myWorkUnit
 *            {job,precursor,fragments[{mz,intensity}],peptides[{id,sequence,mods[{id,mass,residues,num}]}],fixedMods[{id,mass,residues}],fragTol,fargTolUnit}
 * 
 * @returns posts a result string
 *          {job,precursor,peptides[{id,ionsMatched,score,mods[{id,position[]}]}]}
 */
function doModSearch(myWorkUnit) {
	console.warn("doModSearch");
	// console.log(JSON.stringify(myWorkUnit.peptides));
	var allPeptideScores = []; // made up of {peptide id, ionsmatched, score
	// mods:[{ id,position[] }]}
	for (var pep = 0; pep < myWorkUnit.peptides.length; pep++) {
		var currPeptide = myWorkUnit.peptides[pep];
		var resObj = getInitialResObj(currPeptide); // fills out a peptide
		// results object with correct peptide and mod ids
		// {id,ionsMatched,score,mods[{id,position[]}]}

		var currScoreObj;
		if (phosphoModExpected(currPeptide)) {

			// we need to do a search for 2 types of phospho loss as well as
			// modification
			var dum1peptide = getNeutralLossPeptide(currPeptide, PHOSPHO_LOSS);
			var dum2peptide = getNeutralLossPeptide(currPeptide, PHOSPHO2_LOSS);

			var score0 = scorePeptide(currPeptide);
			var score1 = scorePeptide(dum1peptide);
			var score2 = scorePeptide(dum2peptide);

			if (score2.score > score1.score && score2.score > score0.score) {
				currScoreObj = score2;
			} else {
				if (score1.score > score0.score) {
					currScoreObj = score1;
				} else {
					currScoreObj = score0;
				}
			}
		} else {
			currScoreObj = scorePeptide(currPeptide); // currScoreObj is
			// best score for this peptide/mod combo... in form
			// {modPos:[],ionsMatched:0,score:0}
		}

		resObj.score = currScoreObj.score;
		resObj.ionsMatched = currScoreObj.ionsMatched;
		for (var i = 0; i < currScoreObj.modPos.length; i++) { // //currScoreObj.modPos
			// is an array
			// of ModLocs
			var m = currScoreObj.modPos[i].modIndex;
			var p = currScoreObj.modPos[i].possLoc;
			resObj.mods[m].position.push(p); // to resObj.mod[].position;
		}

		// For the moment I have to bulk out the position field to be equal to
		// number of mods (if >=200 then I haven't found/cannot confirm a
		// position)
		for (i = 0; i < resObj.mods.length; i++) {
			if (resObj.mods[i].position.length < currPeptide.mods[i].num) {
				var n = currPeptide.mods[i].num
						- resObj.mods[i].position.length;
				for (var t = 0; t < n; t++) {
					resObj.mods[i].position.push(200 + t);
				}
			}
		}

		// console.log(JSON.stringify(resObj));
		allPeptideScores.push(resObj);
	}
	allPeptideScores.sort(function(a, b) {
		return b.score - a.score;
	});

	var resultObject = {
		job : myWorkUnit.job,
		precursor : myWorkUnit.precursor,
		peptides : []
	}; // what is sent back as retString;

	var PHASE2_MAX_NUMBER_SCORES = 75;
	var outputCount = (myWorkUnit.peptides.length < PHASE2_MAX_NUMBER_SCORES) ? myWorkUnit.peptides.length
			: PHASE2_MAX_NUMBER_SCORES;
	for (var r = 0; r < outputCount; r++) {
		if (allPeptideScores[r].score > 0) {
			resultObject.peptides.push(allPeptideScores[r]);
		}
	}

	// var instring = JSON.stringify(myWorkUnit);
	// console.log("Work Unir: "+instring);
	var retString = JSON.stringify(resultObject);
	console.log("doModSearch return string = " + retString);
	console.log("doModSearch pre = " + resultObject.precursor);

	postMessage(retString);
}

/**
 * 
 * @param currPeptide
 *            peptide object in form {id,sequence,mods[{id,mass,residues,num}]}
 * @returns currScoreObj the best score for this peptide in form
 *          {modPos:[],ionsMatched:0,score:0}
 */
function scorePeptide(currPeptide) {
	var scoreObj = {};
	var currScoreObj = {
		modPos : [],
		ionsMatched : 0,
		score : 0
	};

	// holds the best so far score for this peptide modification )positions
	// held in modpos[{pos,index,mass}]
	var totalModNum = getTotalModNum(currPeptide);

	// Array of ModLoc objects (all possible locations of all modifications
	// reported) ModLoc objects are {possLoc,modIndex,vModMass}]
	var modLocs = getAllModLocs(currPeptide);

	// just check..
	if (modLocs.length < totalModNum) {
		// we are looking for more mods than possible with this residue because
		// we cheat and place STY when ANDREW's data has full complement.
		console.log("BEWARE DATA ERROR");
		// even worse now as we look for 21 twice
		totalModNum = modLocs.length;
	}

	// 24/5/17 Boss wants me to try with the expected number of mods
	var subIonsets = getIonsetForAllMods(currPeptide, modLocs, totalModNum);

	// an array of ionsets to score containing num number of mods
	// console.log("ionsets="+JSON.stringify(subIonsets));
	for (var s = 0; s < subIonsets.length; s++) {
		subIonsets[s] = matchSpectraWithIonSet(g_myWorkUnit.fragments,
				subIonsets[s]); // for each ionset we log matches with ms2
		// fragments
		scoreObj = scoreIonset(subIonsets[s]); // score the ionset (matched
		// ions and ion ladders)
		/*
		 * if (currPeptide.id === 1021615){
		 * console.log(JSON.stringify(subIonsets[s])); console.log("Id:
		 * "+currPeptide.id+" "+currPeptide.sequence+"
		 * modPos:"+JSON.stringify(subIonsets[s].modPos)+" gives
		 * "+JSON.stringify(scoreObj)); var matchString = ""; var modString =
		 * ""; for (var t = 0; t < subIonsets[s].bIons.length; t++){ if (num >
		 * 0){ modString +=(subIonsets[s].modPos[0].possLoc==(t+1))?"M":"0"; }
		 * matchString +=(subIonsets[s].bIons[t].match)?"X":"0"; } matchString
		 * +="\n"; for (var v = subIonsets[s].yIons.length-1; v >=0; v--){
		 * matchString +=(subIonsets[s].yIons[v].match)?"X":"0"; }
		 * modString+="\n"; matchString +="\n";
		 * console.log(modString+matchString+"-------------"); }
		 */

		if (scoreObj.score >= currScoreObj.score) {
			// CHANGE FROM > TO >= and matched mascot better!!!!!
			currScoreObj.score = scoreObj.score;
			currScoreObj.modPos = subIonsets[s].modPos.slice();
			// place modPos structure [ModLoc] is kept with score
			currScoreObj.ionsMatched = scoreObj.ionsMatched;
		}
	}

	if (currScoreObj.ionsMatched > 5) {
		console.log("Id: " + currPeptide.id + " " + currPeptide.sequence
				+ " mod: " + currPeptide.mods[0].id + " at Pos: "
				+ JSON.stringify(currScoreObj.modPos[0].possLoc) + " score: "
				+ currScoreObj.score);
		console.log("------------");
	}

	return currScoreObj;
}

function getNeutralLossPeptide(peptide, loss) {
	var dumPeptide = {
		id : peptide.id,
		sequence : peptide.sequence,
		mods : [ {
			id : 0,
			mass : 0,
			residues : "",
			num : 1
		} ]
	};

	for (var i = 0; i < peptide.mods.length; i++) {
		dumPeptide.mods[i].id = peptide.mods[i].id;
		dumPeptide.mods[i].mass = peptide.mods[i].mass;
		dumPeptide.mods[i].residues = peptide.mods[i].residues;
		dumPeptide.mods[i].num = peptide.mods[i].num;

		if (dumPeptide.mods[i].id === PHOSPHORYLATION) {
			dumPeptide.mods[i].mass = PHOSPHO_MASS - loss;
		}
	}

	return dumPeptide;
}
/**
 * 
 * @param peptide
 *            Peptide object in form {id,sequence,mods[{id,mass,residues,num}]}
 * @returns true if a PHOSPHORYLATION modification is noted in mods array
 */
function phosphoModExpected(peptide) {
	for (var i = 0; i < peptide.mods.length; i++) {
		if (peptide.mods[i].id === PHOSPHORYLATION) {
			return true;
		}
	}

	return false;
}

if (typeof Math.log10 == "undefined") {
	Object.prototype.log10 = function(n) {
		return Math.log(n) / Math.log(10);
	};
}

function doThirdPhaseSearch(myWorkUnit) {
	var dumPhospho1 = {
		id : 0,
		sequence : "",
		mod : [ {
			id : 21,
			mass : PHOSPHO_MASS - PHOSPHO_LOSS,
			residues : "STY",
			num : 1
		} ]
	};

	var allPeptideScores = []; // made up of {peptide id, ionsmatched, score
	// mods:[{ id,position[] }]}
	for (var pep = 0; pep < myWorkUnit.peptides.length; pep++) {
		var scoreObj = {};
		var currScoreObj = {
			modPos : [],
			ionsMatched : 0,
			score : 0
		}; // holds the best so far score for this peptide modification
		// )positions held in modpos[{pos,index,mass}]

		var currPeptide = myWorkUnit.peptides[pep];

		var totalModNum = getTotalModNum(currPeptide);
		var modLocs = getAllModLocs(currPeptide); // Array of ModLoc objects
		// (all possible locations
		// of all modifications
		// reported) ModLoc objects
		// are
		// {possLoc,modIndex,vModMass}]
		var resObj = getInitialResObj(currPeptide); // flushes out resObj arrays
		// to correct size

		// just check..
		if (modLocs.length < totalModNum) { // we are looking for more mods than
			// possible with this residue
			// because we cheat and place STY
			// when ANDREW's data has full
			// complement.
			// console.log("BEWARE DATA ERROR")
			totalModNum = modLocs.length; // even worse now as we look for 21
			// twice
		}
		// 24/5/17 Boss wants me to try with the expected number of mods...
		// hence the +totalModNum.
		// I used to check for zero mods being there to totalModNum as I need to
		// id the underlying peptide.
		for (var num = 0 + totalModNum; num <= totalModNum; num++) {
			// as some mods fall off we are going to check for 0 mods to num
			// mods
			var subIonsets = getIonsetForAllMods(currPeptide, modLocs, num); // an
			// array of ionsets to score containing num number of mods ionsets
			// look good.
			// console.log("ionset0="+JSON.stringify(subIonsets[0]));

			for (var s = 0; s < subIonsets.length; s++) {
				subIonsets[s] = matchSpectraWithIonSet(myWorkUnit.fragments,
						subIonsets[s]); // for each ionset we log matches with
				// ms2 fragments
				scoreObj = scoreIonset(subIonsets[s]);
				// score the ionset (matched ions and ion ladders)

				if (scoreObj.score >= currScoreObj.score) {
					// CHANGE FROM > TO >= and matched mascot better!!!!!
					currScoreObj.score = scoreObj.score;
					currScoreObj.modPos = subIonsets[s].modPos.slice();
					// place modPos structure [ModLoc] is kept with score
					currScoreObj.ionsMatched = scoreObj.ionsMatched;
				}
			}
		}

		// currScoreObj = best score, ion-match and modPos of this best score
		resObj.score = currScoreObj.score;
		resObj.ionsMatched = currScoreObj.ionsMatched;
		// sort out resobj mod positions
		for (var i = 0; i < currScoreObj.modPos.length; i++) {
			// currScoreObj.modPos is an array of ModLocs
			var m = currScoreObj.modPos[i].modIndex;
			var p = currScoreObj.modPos[i].possLoc;
			resObj.mods[m].position.push(p); // to resObj.mod[].position;
		}

		// For the moment I have to bulk out the position field to be equal to
		// number of mods (if >=200 then I haven't found/cannot confirm a
		// position)
		for (var index = 0; index < resObj.mods.length; index++) {
			if (resObj.mods[index].position.length < currPeptide.mods[index].num) {
				var n = currPeptide.mods[index].num
						- resObj.mods[index].position.length;
				for (var t = 0; t < n; t++) {
					resObj.mods[index].position.push(200 + t);
				}
			}
		}
		allPeptideScores.push(resObj);
	}

	allPeptideScores.sort(function(a, b) {
		return b.score - a.score;
	});

	var resultObject = {
		job : myWorkUnit.job,
		precursor : myWorkUnit.precursor,
		peptides : []
	}; // what is sent back as retString;

	var PHASE2_MAX_NUMBER_SCORES = 75;
	var outputCount = (myWorkUnit.peptides.length < PHASE2_MAX_NUMBER_SCORES) ? myWorkUnit.peptides.length
			: PHASE2_MAX_NUMBER_SCORES;
	for (var r = 0; r < outputCount; r++) {
		// if (allPeptideScores[r].score > 0){
		resultObject.peptides.push(allPeptideScores[r]);
		// }
	}

	if (typeof WorkerGlobalScope !== 'undefined'
			&& self instanceof WorkerGlobalScope) {
		postMessage(resultObject);
		return;
	}

	sendResult(resultObject);
}

function getTotalModNum(myPeptide) {
	var mnum = 0;
	for (var mod = 0; mod < myPeptide.mods.length; mod++) {
		// we may have more than one possible modification per peptide
		mnum += myPeptide.mods[mod].num;
	}

	return mnum;
}

function getInitialResObj(myPeptide) {
	var resObj = {
		id : myPeptide.id,
		ionsMatched : 0,
		score : 0,
		mods : []
	};

	if (myPeptide.mods === undefined) {
		return resObj;
	}

	for (var mod = 0; mod < myPeptide.mods.length; mod++) {
		var modp = {
			id : myPeptide.mods[mod].id,
			position : []
		};

		resObj.mods.push(modp);
		// flush out resObj.mods[]
	}

	return resObj;
}

function getAllModLocs(myPeptide) {
	var mlocs = [];
	if (myPeptide.mods === undefined)
		return mlocs;
	for (var mod = 0; mod < myPeptide.mods.length; mod++) {
		// we may have more than one possible modification per peptide FUDGE
		if (myPeptide.mods[mod].id === PHOSPHORYLATION) {
			myPeptide.mods[mod].residues = "STY";
			// myPeptide.mods[mod].mass -= PHOSPHO_LOSS;
			// lets just do phospho loss!
		}
		// END FUDGE
		var possLocs = createIndexofPossibleLocations3(myPeptide.sequence,
				myPeptide.mods[mod].residues);

		for (var i = 0; i < possLocs.length; i++) {
			var modObj = new ModLoc(possLocs[i], mod, myPeptide.mods[mod].mass);
			mlocs.push(modObj);
		}
	}
	mlocs.sort(function(a, b) {
		return a.possLoc - b.possLoc;
	}); // why do we sort it? to shuffle the mods?

	return mlocs;
}

function getIonsetForAllMods(myPeptide, modlocs, num) {
	// array ionsets to return
	var ionSetArray = [];

	// combination array when dealing with more than one number of modifications
	// of one type
	var combArray = [];

	// /var modMass = myMod.mass //shorthand;

	var MUCH_TOO_MUCH = 10000; // Arbitrary value at which point it becomes to
	// large causing memory issues and to long to
	// process (Server will assume dropped
	// connection)

	// search for none modified peptide - may well be here (Peaks PTM) and gives
	// a basal score for peptide.
	if (num === 0) {
		ionSetArray.push(getIonsFromArray3(myPeptide, [], 0)); // single ionset
		// returned

		return ionSetArray;
	}

	// search for occurrence of a single modification at any of the possible
	// possible locations
	// the others may have fallen off. Finding a good ion match can reduce size
	// of future search
	if (num === 1) {
		for (var ml = 0; ml < modlocs.length; ml++) {
			ionSetArray.push(getIonsFromArray3(myPeptide, [ modlocs[ml] ]));
		}

		return ionSetArray;
	}

	if (modlocs.length < num) {
		return ionSetArray;
	}
	// do a calculation first to see if do-able combinations rapidly get out of
	// hand.
	// number of combinations = n!/(k!(n-k)!)
	if (getCombinations(modlocs.length, num) > MUCH_TOO_MUCH) {
		console.log("Bailing: Array too large for me.");
		return ionSetArray;
	}

	combArray = createArrayPossibleCombinations3(modlocs, modlocs.length, num); // combarray
	// now
	// in
	// form
	// of
	// [modlocs]

	for (var c = 0; c < combArray.length; c++) {
		var mlocs = combArray[c]; // mlocs = confirmed mods + combination to
		// try
		// console.log("mlocs:"+JSON.stringify(mlocs)+" len:"+modlocs.length+"
		// num:"+num+"conpos:"+conpos.length);
		var ions = getIonsFromArray3(myPeptide, mlocs);

		ionSetArray.push(ions);
	}

	return ionSetArray;
}

function getIonsFromArray3(myPeptide, mlocs) // looking for more than one
// modification at a time
{
	var ionSet = new Ionset(); // {modPos:[],bIons:[],yIons:[]}; //bIons =
	// [{mass, match, intensity}]
	var acid = ''; // amino acid character
	var sequence = myPeptide.sequence; // 'ABCEF'

	ionSet.modPos = mlocs.slice(); // initialise modpos with the array od
	// {possloc,index,mass}

	var cumulativeMass = 1.007276; // cumulative mass for b ions start with
	// Water - OH which is H
	cumulativeMass += checkforFixedPTM('['); // is a fixed mod an amine
	// terminus?
	var ionObj;
	var loopIndex;
	for (var b = 0; b < (sequence.length); b++) {
		ionObj = new Ion(); // {mass:0,match:0,intensity:0,deltaM:0,modFlag:false};
		acid = sequence.charAt(b);
		cumulativeMass += g_AAmass[acid] + checkforFixedPTM(acid);

		for (loopIndex = 0; loopIndex < mlocs.length; loopIndex++) {
			if (b === (mlocs[loopIndex].possLoc - 1)) {
				cumulativeMass += mlocs[loopIndex].vModMass;
				ionObj.modFlag = true;
			}
		}

		ionObj.mass = cumulativeMass;
		ionSet.bIons.push(ionObj);
	}

	cumulativeMass = 18.010565 + 1.007276; // start with H + water
	cumulativeMass += checkforFixedPTM(']'); // is a fixed mod a carboxyl
	// end?

	for (var y = sequence.length - 1; y >= 0; y--) {
		ionObj = new Ion(); // {mass:0,match:0,intensity:0,deltaM:0,modFlag:false};
		acid = sequence.charAt(y);
		cumulativeMass += g_AAmass[acid] + checkforFixedPTM(acid);

		for (loopIndex = 0; loopIndex < mlocs.length; loopIndex++) {
			if (y === (mlocs[loopIndex].possLoc - 1)) {
				cumulativeMass += mlocs[loopIndex].vModMass;
				ionObj.modFlag = true;
			}
		}
		ionObj.mass = cumulativeMass;
		ionSet.yIons.push(ionObj);
	}

	return ionSet;
}

// returns integer array of all possible locations of mod +1
// NB position = sequence position + 1 ... 0 = '[' 1 to len = = to len-1 and
// length+1 = ] ]
// input sequence (eg "PEPTIDE") and residues (eg "CKV");
function createIndexofPossibleLocations3(sequence, residues) {
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
			pos = residues.length + 1;
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
}

// -------------------------------------------------------
// -------------------------------------------------------
// adds unique parts of second to first until maxlen is reached
function uniqueArrayConcat(first, second, maxlen) {
	if (first.length >= maxlen) {
		return first;
	}

	for (var s = 0; s < second.length; s++) { // check new locations held in
		// contemp[] to see if already
		// in conpos
		var val = second[s];
		var valmatch = false;
		for (var f = 0; f < first.length; f++) {
			if (first[f] === val) {
				valmatch = true; // already in conpos so no need to add this
			}
		}

		if (!valmatch) {
			if (first.length < maxlen) {
				first.push(val); // new one found in this ions set so add to
				// end;
			}
		}
	}
	return first;
}

// -----------------------------------------------------------------------------
// Combination & Factorial code

function createArrayPossibleCombinations3(locArray, n, k) {
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
}

function getCombinations(n, k) {
	var nfac = getFactorial(n);
	var kfac = getFactorial(k);
	var nsubkfac = getFactorial(n - k);

	return nfac / (kfac * nsubkfac);
}

function getFactorial(f) {
	var r = 1;
	for (; f > 0; r *= f, f--) {
	}

	return r;
}