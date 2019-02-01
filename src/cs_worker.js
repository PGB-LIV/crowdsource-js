"use strict";

/**
 * Copyright 2018 University of Liverpool
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

/**
 * Web Worker
 */
function doSearch(data) {
	g_myWorkUnit = data;

	// includes variable modifications
	doThirdPhaseSearch(g_myWorkUnit);
}

this.onmessage = function(event) {
	doSearch(event.data);
};
