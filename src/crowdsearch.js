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
var BLOB = true; // true if cross domian required

/**
 * The PHP that handles data requests and accepts results (will need to go to
 * master)
 */
var PROVIDER_PHP = "http://pgb.liv.ac.uk/~andrew/crowdsource-server/src/public_html/job.php";

/**
 * WebWorker instance
 */
var myWorker;
var myWorkUnit;

/**
 * Security constraints stop us building a cross domain worker. We need to add
 * it to our instance via a Blob URL
 */
var BuildWorker = function(foo) {
	var str = foo.toString().match(
			/^\s*function\s*\(\s*\)\s*\{(([\s\S](?!\}$))*[\s\S])/)[1];

	return new Worker(window.URL.createObjectURL(new Blob([ str ], {
		type : 'text/javascript'
	})));
};

if (typeof (Worker) == "undefined") {
	throw new Error("WebWorker unsupported");
}

initialiseWorker();
console.log("myWorker = " + myWorker);
requestWorkUnit();
console.log("requested first Work unit");

function initialiseWorker() {
	if (BLOB) {
		myWorker = BuildWorker(function() {
			var WORKER_LOCATION = "http://pgb.liv.ac.uk/~johnheap/crowdsource-server/src/public_html/javascript/";
			// will need to go to master
			importScripts(WORKER_LOCATION + "cs_worker.js");
			importScripts(WORKER_LOCATION + "thirdphase_worker.js");

		});
	} else {
		myWorker = new Worker(
				"http://pgb.liv.ac.uk/~johnheap/crowdsource-server/src/public_html/javascript/no_blob_worker.js");
	}
}

function receiveWorkUnit(json) {
	console.log("received a workunit");
	myWorkUnit = JSON.stringify(json);
	console.log("received a workunit: " + myWorkUnit);
	// on receiving workUnit it will get to work
	myWorker.postMessage(myWorkUnit);
}

// worker communicates with the main js via JSON strings
myWorker.onmessage = function(e) {
	var workerResponse = JSON.parse(e.data);
	if (typeof (workerResponse.job) !== "undefined") {
		sendResult(JSON.stringify(workerResponse));
	}
};

// terminate session event. Also occurs on Refresh.
$(window).on("beforeunload", function() {
	if (typeof (w) !== "undefined") {
		myWorker.terminate(); // terminate the worker...
		sendTerminating();
	}
});

/**
 * Server Communication
 */

function requestWorkUnit() {
	$.getScript(PROVIDER_PHP + "?r=workunit");
}

function sendResult(resultString) {
	$.getScript(PROVIDER_PHP + "?r=result&result=" + resultString);
}

function sendTerminating() {
	$.getScript(PROVIDER_PHP + "?r=terminate");
}

/**
 * this function is the P of the JSONP response from server eg
 * "parseResult(Object)" P the response server
 */
function parseResult(json) {
	if (typeof (json.job) !== "undefined") {
		// If job is defined, work unit has been sent
		receiveWorkUnit(json);
		return;
	}

	switch (json.type) {
	case "nomore":
		console.log("All done.");
		break;
	case "confirmation":
		requestWorkUnit();
		break;
	default:
		break;
	}
}