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

/**
 * GLOBAL {job:0, precursor:0, fragments:[{mz:n, intensity:n}...],
 * peptides:[{id:1, structure:"ASDFFS"}...]};
 */

var g_myWorkUnit;
var g_AAmass = {
	A : 71.037114,
	R : 156.101111,
	N : 114.042927,
	D : 115.026943,
	C : 103.009185,
	E : 129.042593,
	Q : 128.058578,
	G : 57.021464,
	H : 137.058912,
	I : 113.084064,
	L : 113.084064,
	K : 128.094963,
	M : 131.040485,
	F : 147.068414,
	P : 97.052764,
	S : 87.032028,
	T : 101.047679,
	U : 150.95363,
	W : 186.079313,
	Y : 163.063329,
	V : 99.068414,
	X : 0
};

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
	this.match = false;
	this.intensity = 0;
	this.deltaM = 0;
	this.modFlag = false;
}

/**
 * used in second/third phase
 */
function ModLoc(ploc, modindex, modmass) {
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
}

var ION_MATCH_SCALE = 20;
var ION_DELTA_M_SCALE = 0;
var ION_LADDER_SCALE = 30;
var ION_INTENSITY_SCALE = 2;

function scoreIonset(ionSet) {
	var scoreObj = {
		ionsMatched : 0,
		score : 0
	};

	var bcount = 0;
	var bintensity = 1;
	var bscore = 0;
	var ycount = 0;
	var yscore = 0;
	var yintensity = 1; /* as log 0 = null */

	for (var b = 0; b < ionSet.bIons.length; b++) {
		if (ionSet.bIons[b].match) {
			bcount++;
			bintensity += ionSet.bIons[b].intensity;
		}
	}

	// 0 - 25% for example.
	bscore = (ION_MATCH_SCALE * bcount) / ionSet.bIons.length;
	bscore += (ION_INTENSITY_SCALE * Math.log10(bintensity) / ionSet.bIons.length);

	for (var y = 0; y < ionSet.yIons.length; y++) {
		if (ionSet.yIons[y].match) {
			ycount++;
			yintensity += ionSet.yIons[y].intensity;
		}
	}

	// 0 -25 for example
	yscore = (ION_MATCH_SCALE * ycount) / ionSet.yIons.length;
	yscore += (ION_INTENSITY_SCALE * Math.log10(yintensity) / ionSet.yIons.length);

	// the idea being that true discoveries will have more intact ion ladders
	var consecResult = check_consecutive_ions(ionSet);
	bscore += (ION_LADDER_SCALE * consecResult.bconseq) / ionSet.bIons.length;
	yscore += (ION_LADDER_SCALE * consecResult.yconseq) / ionSet.yIons.length;

	scoreObj.score = yscore + bscore;
	/*
	 * we have modified the score in all cases by dividing by the peptide length
	 * //so 4 counts of 8 res peptide is better than 4 counts of 16 res
	 * peptide... //but 8 counts of 16 res peptide is notionally better than 4
	 * counts of 8 res peptide. Can we also modifiy the score by a samll factor
	 * based on length. //try
	 */

	scoreObj.score = scoreObj.score * (1 + (ionSet.bIons.length / 120));

	scoreObj.score = scoreObj.score.toFixed(2);
	scoreObj.score *= 1;
	scoreObj.ionsMatched = bcount + ycount;
	return scoreObj;
}

/*
 * Differentiation between real and decoy data is enhanced by scoring
 * consecutive ions //This function doesn't take into account the length of the
 * peptide - this is adjusted for by the calling function scoreIonset(ionSet)
 */
function check_consecutive_ions(ionSet) {
	var MININUM_CONSECUTIVE_COUNT = 3; // ignore if less than this
	// number of consecutive ions
	var bconseq = 0;
	var yconseq = 0;
	var bcnt = 0;
	var ycnt = 0;
	var consecResult = {
		bconseq : 0,
		yconseq : 0
	};

	for (var b = 0; b < ionSet.bIons.length; b++) {
		if (ionSet.bIons[b].match) {
			bcnt++;
			if (bcnt > bconseq) {
				bconseq = bcnt;
			}
		} else {
			bcnt = 0;
		}
	}

	for (var y = 0; y < ionSet.yIons.length; y++) {
		if (ionSet.yIons[y].match) {
			ycnt++;
			if (ycnt > yconseq) {
				yconseq = ycnt;
			}
		} else {
			ycnt = 0;
		}
	}

	if (bconseq < MININUM_CONSECUTIVE_COUNT) {
		bconseq = 0;
	}

	if (yconseq < MININUM_CONSECUTIVE_COUNT) {
		yconseq = 0;
	}

	consecResult.bconseq = bconseq;
	consecResult.yconseq = yconseq;

	return consecResult;
}

// common
function matchSpectraWithIonSet(spectra, ionSet, checkloss) {
	if (checkloss === undefined) {
		checkloss = true;
	}

	var resultObject = {
		deltaM : 0,
		intensity : 0
	}; // modifying
	for (var b = 0; b < ionSet.bIons.length; b++) {
		resultObject = massFoundInSpectra(spectra, ionSet.bIons[b].mass,
				checkloss);
		ionSet.bIons[b].deltaM = resultObject.deltaM;
		ionSet.bIons[b].intensity = resultObject.intensity;
		if (ionSet.bIons[b].intensity !== 0) {
			ionSet.bIons[b].match = true; // could remove this and
			// use non zero
			// intensity as match
			// boolean
		}
	}

	for (var y = 0; y < ionSet.yIons.length; y++) {
		resultObject = massFoundInSpectra(spectra, ionSet.yIons[y].mass,
				checkloss);
		ionSet.yIons[y].intensity = resultObject.intensity;
		ionSet.yIons[y].deltaM = resultObject.deltaM;
		if (ionSet.yIons[y].intensity !== 0) {
			ionSet.yIons[y].match = true;
		}
	}

	return ionSet; // we've modified it I know don't need to return it
	// as modified by function but it makes it explicit.
}

// Is ion mass (mass) found in ms2 spectrum (spectra) returns zero if no match
// else first match intensity. also store deltaM (how close to precision window
// we are) common
function massFoundInSpectra(spectra, mass, checkloss) {
	checkloss = false; // AMMONIA and WATER loss seem to be red
	// herrings and worsen results.
	var WATER = 18.010565;
	var AMMONIA = 17.026549;
	var resultObject = {
		deltaM : 0,
		intensity : 0
	};
	var precision = (g_myWorkUnit.fragTolUnit === "ppm") ? mass * 0.000001
			* g_myWorkUnit.fragTol : g_myWorkUnit.fragTol;
	var deltaM = 0;

	for (var m = 0; m < spectra.length; m++) {
		// NB file format we are using mz is actually neutral charge. (Must
		// check if change file type)
		deltaM = Math.abs(mass - spectra[m].mz);
		if (deltaM <= precision) {
			if (spectra[m].intensity > resultObject.intensity) {
				resultObject.deltaM = 1 - (deltaM / precision);
				resultObject.intensity = spectra[m].intensity;
			}
		}

		if (checkloss) {
			deltaM = Math.abs((mass - WATER) - spectra[m].mz);
			if (deltaM <= precision) {
				if (spectra[m].intensity > resultObject.intensity) {
					resultObject.deltaM = 1 - (deltaM / precision);
					resultObject.intensity = spectra[m].intensity;
				}
			}

			deltaM = Math.abs((mass - AMMONIA) - spectra[m].mz);
			if (deltaM <= precision) {
				if (spectra[m].intensity > resultObject.intensity) {
					resultObject.deltaM = 1 - (deltaM / precision);
					resultObject.intensity = spectra[m].intensity;
				}
			}
		}
	}
	return resultObject;
}

function checkforFixedPTM(res) {
	var masstoadd = 0;
	for (var m = 0; m < g_myWorkUnit.fixedMods.length; m++) {
		if (g_myWorkUnit.fixedMods[m].residues.length > 1) { // just to
			// check...
			console.log("Fixed Mod needs looking at "
					+ g_myWorkUnit.fixedMods[m].residues);
		}
		if (g_myWorkUnit.fixedMods[m].residues === res) {
			masstoadd += g_myWorkUnit.fixedMods[m].mass;
			// console.log("Adding "+masstoadd+" to
			// "+g_myWorkUnit.mods[m]['loc']);
		}
	}

	return masstoadd;
}

/**
 * Web Worker
 */
function doSearch(data) {
	g_myWorkUnit = data;

	// change the defaults if phase2/3 ion matches are more important.
	ION_MATCH_SCALE = 30;
	ION_LADDER_SCALE = 20;

	// includes variable modifications
	doThirdPhaseSearch(g_myWorkUnit);
}

this.onmessage = function(event) {
	doSearch(event.data);
};
