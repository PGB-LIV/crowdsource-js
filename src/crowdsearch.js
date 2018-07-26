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
var BLOB = false; // true if corss domian required

/**
 * The PHP that handles data requests and accepts results (will need to go to
 * master)
 */
var PROVIDER_PHP = "http://pgb.liv.ac.uk/~andrew/crowdsource-server/src/public_html/job.php"; // 

var myWorker; // WebWorker instance
var myWorkUnit;

// Security constraints stop us building a cross domain worker. We need to add
// it to our instance via a Blob URL
var BuildWorker = function(foo) {
	var str = foo.toString().match(
			/^\s*function\s*\(\s*\)\s*\{(([\s\S](?!\}$))*[\s\S])/)[1];
	return new Worker(window.URL.createObjectURL(new Blob([ str ], {
		type : 'text/javascript'
	})));
};

if (typeof (Worker) == "undefined") {
	console.log("NO WEB WORKER SUPPORT!!!!!!!!!!!!!!!!!!!!!!");
}
initialiseWorker();
console.log("myWorker = " + myWorker);
requestWorkUnit();
console.log("requested first Work unit");

function initialiseWorker() {
	if (BLOB) {
		myWorker = BuildWorker(function() {
			var WORKER_LOCATION = "http://pgb.liv.ac.uk/~johnheap/crowdsource-server/src/public_html/javascript/"; // will
			// need
			// to
			// go
			// to
			// master
			var WORKER1_NAME = "cs_worker.js"; // name of worker javascript
			var WORKER2_NAME = "firstphase_worker.js";
			var WORKER3_NAME = "secondphase_worker.js";

			importScripts(WORKER_LOCATION + WORKER1_NAME);
			importScripts(WORKER_LOCATION + WORKER2_NAME);
			// importScripts(WORKER_LOCATION+WORKER3_NAME);
			importScripts(WORKER_LOCATION + "thirdphase_worker.js");

		});
	} else {
		myWorker = new Worker(
				"http://pgb.liv.ac.uk/~johnheap/crowdsource-server/src/public_html/javascript/"
						+ 'no_blob_worker.js');
	}
}

function requestWorkUnit() {
	$.getScript(PROVIDER_PHP + "?r=workunit");
}

function receiveWorkUnit(json) {
	console.log("received a workunit");
	myWorkUnit = JSON.stringify(json);
	console.log("received a workunit: " + myWorkUnit);
	myWorker.postMessage(myWorkUnit); // on receiving workUnit it will get to
	// work
}

function sendResult(resultString) {
	$.getScript(PROVIDER_PHP + "?r=result&result=" + resultString);
}

function sendTerminating() {
	$.getScript(PROVIDER_PHP + "?r=terminate");
}

// this function is the P of the JSONP response from server eg
// "parseResult(Object)"
function parseResult(json) {
	if (typeof (json.job) !== "undefined") {
		receiveWorkUnit(json); // there is no type in workunits no more. so
		// check if job is defined
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

myWorker.onmessage = function(e) // worker communicates with the main js via
// JSON strings
{
	var workerResponse = JSON.parse(e.data);
	if (typeof (workerResponse.job) !== "undefined") {
		sendResult(JSON.stringify(workerResponse));
	}
}

// terminate session event. Also occurs on Refresh.
$(window).on("beforeunload", function() {
	if (typeof (w) !== "undefined") {
		myWorker.terminate(); // terminate the worker...
		sendTerminating();
	}
});
