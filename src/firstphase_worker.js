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

// This file holds first phase specific functions.

/*
 * NOW AFTER REFACTORING myWorkUnit =
 * {"job":1,"precursor":1,"fragments":[{"mz":102.055527,"intensity":12677.620117},..],
 * "peptides":[{"id":44982,"sequence":"NHFNANCSFWNYSER",varmod:{id,
 * mass,"residue sequence"}},...]
 * "fixedMods":[{"mass":57.021464,"residue":"C"}],"fragTol":0.01,"fragTolUnit":"da"}
 * "varmod" is optional and occurs in pass 2 type workunits
 */

// First phase search only deals with fixed mods - so only one set of ion
// fragments to check
function doFirstPhaseSearch(myWorkUnit) {
	var allPeptideScores = []; // array of scoreObjects {id:0, score:0} one per
								// peptide matched
	var resultObject = {
		job : myWorkUnit.job,
		precursor : myWorkUnit.precursor,
		peptides : []
	}; // return object
	var scoreObj = {
		id : 0,
		ionsMatched : 0,
		score : 0
	};
	var numRes = MAX_NUMBER_SCORES_RETURNED; // number of results to return

	for (var p = 0; p < myWorkUnit.peptides.length; p++) {
		var ionset = fixedOnlyFragment(myWorkUnit.peptides[p]); // ionset =
																// {modpos[],
																// bIons[],yIons[]}
																// //we don't
																// need modpos
																// in this phase
		ionset = matchSpectraWithIonSet(myWorkUnit.fragments, ionset); // fills
																		// out
																		// match
																		// and
																		// intensity
																		// in
																		// bIons[]
																		// and
																		// yIons[];
		scoreObj = scoreIonset(ionset); // fills out score and ionsMatched in
										// {ionsMatched:0, score: 0};
		scoreObj.id = myWorkUnit.peptides[p].id; // add an id
		allPeptideScores.push(scoreObj);
	}
	allPeptideScores.sort(function(a, b) {
		return b.score - a.score
	});

	if (myWorkUnit.peptides.length < numRes) {
		numRes = myWorkUnit.peptides.length;
	}
	for (var r = 0; r < numRes; r++) {
		if (allPeptideScores[r].score > 0) {
			resultObject.peptides.push(allPeptideScores[r])
		}
	}
	var retString = JSON.stringify(resultObject);
	// console.log(retString);
	console.log("job = " + resultObject.job + ", pre = "
			+ resultObject.precursor + ", peps tested = "
			+ myWorkUnit.peptides.length + ", peps returned = "
			+ resultObject.peptides.length);
	postMessage(retString);
}

// Phase one fragmentation. - only fixed modifications are taken into account.
// returns a single phase one ionSet object
// {bIons:[{mass,match,intensity,deltaM},..],yIons:[{mass,match,intensity,deltaM},..]}
function fixedOnlyFragment(myPeptide) {
	var ionSet = new Ionset(); // {bIons:[],yIons:[]}; //bIons = [{mass, match,
								// intensity}]
	var acid = ''; // amino acid character
	// var acidMass = 0;
	// var fixModMass = 0; //mass to add to this residue due to a fixed mod
	var sequence = myPeptide.sequence; // 'ABCEF'

	var cumulativeMass = 1.007276; // cumulative mass for b ions start with
									// Water - OH which is H
	cumulativeMass += checkforFixedPTM('['); // is a fixed mod an amine
												// terminus?
	for (var b = 0; b < (sequence.length); b++) {
		var ionObj = new Ion(); // {mass:0,match:0,intensity:0,deltaM:0};
		acid = sequence.charAt(b);
		// fixModMass = checkforFixedPTM(acid);
		// acidMass = g_AAmass[acid]+fixModMass;
		// cumulativeMass+=acidMass;
		cumulativeMass += g_AAmass[acid] + checkforFixedPTM(acid);
		ionObj.mass = cumulativeMass;
		ionSet.bIons.push(ionObj);
	}

	cumulativeMass = 18.010565 + 1.007276; // start with H + water
	cumulativeMass += checkforFixedPTM(']'); // is a fixed mod a carboxyl
												// end?
	for (var y = sequence.length - 1; y >= 0; y--) {
		var ionObj = new Ion(); // {mass:0,match:0,intensity:0,deltaM:0};
		acid = sequence.charAt(y);
		// fixModMass = checkforFixedPTM(acid);
		// acidMass = g_AAmass[acid]+fixModMass;
		// cumulativeMass+=acidMass;
		cumulativeMass += g_AAmass[acid] + checkforFixedPTM(acid);
		ionObj.mass = cumulativeMass;
		ionSet.yIons.push(ionObj);
	}
	return ionSet;
}
