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
 * The PHP that handles data requests and accepts results (will need to go to
 * master)
 */
var MASTER_URL = "http://pgb.liv.ac.uk/~andrew/crowdsource-server/src/public_html/job.php";
var SCRIPT_URL = "http://pgb.liv.ac.uk/~andrew/crowdsource-js/src/";
var CALLBACK = "parseResult";
var MAX_JOB_COUNT = 10;
var JOB_COUNT = 0;

/**
 * WebWorker instance
 */
var myWorker;

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

function initialiseWorker() {
	myWorker = BuildWorker(function() {
		// will need to go to master
		importScripts(SCRIPT_URL + "cs_worker.js");
		importScripts(SCRIPT_URL + "thirdphase_worker.js");
	});

	// worker communicates with the main js via JSON strings
	myWorker.onmessage = function(e) {
		sendResult(e.data);
	};
}

function receiveWorkUnit(json) {
	if (myWorker === undefined) {
		doSearch(json);
		return;
	}

	myWorker.postMessage(json);
}

// terminate session event. Also occurs on Refresh.
$(window).on("beforeunload", function() {
	if (typeof (w) !== "undefined") {
		myWorker.terminate(); // terminate the worker...
		sendTerminating();
	}
});

if (typeof (Worker) !== "undefined") {
	$.ajax({
		url : SCRIPT_URL + 'cs_worker.js',
		dataType : 'script',
		async : false
	});

	$.ajax({
		url : SCRIPT_URL + 'thirdphase_worker.js',
		dataType : 'script',
		async : false
	});

	console.info("WebWorker not available");
} else {
	initialiseWorker();
}

requestWorkUnit();

/**
 * Server Communication
 */

function requestWorkUnit() {
	if (JOB_COUNT > MAX_JOB_COUNT) {
		return;
	}

	JOB_COUNT++;
	$.getScript(MASTER_URL + "?r=workunit&callback=" + CALLBACK);
}

function sendResult(resultObject) {
	var resultString = JSON.stringify(resultObject);

	$.getScript(MASTER_URL + "?r=result&result=" + resultString + "&callback="
			+ CALLBACK);
}

function sendTerminating() {
	$.getScript(MASTER_URL + "?r=terminate&callback=" + CALLBACK);
}

/**
 * this function is the P of the JSONP response from server eg
 * "parseResult(Object)" P the response server
 */
function parseResult(json) {
	if (typeof (json.job) !== "undefined") {
		// If job is defined, work unit has been sent

		console.log("Received work unit" + JOB_COUNT);
		receiveWorkUnit(json);
		return;
	}

	switch (json.type) {
	case "confirmation":
		console.log("Requesting work unit");
		requestWorkUnit();
		break;
	default:
	case "nomore":
		break;
	}
}
