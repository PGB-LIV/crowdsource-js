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

function DraculaClient(callBackInstance) {
	"use strict";

	var REQUEST_URI = "http://138.253.218.170:1260/work";
	var WORKER_URI = "http://pgb.liv.ac.uk/~andrew/crowdsource-js/src/MsSearch.js";

	var draculaInstance = this;

	this.callBack = callBackInstance + ".parseResult";

	this.renfield;

	this.jobCount = 0;

	/**
	 * this function is the P of the JSONP response from server eg
	 * "parseResult(Object)" P the response server
	 */
	this.parseResult = function(json) {
		if (typeof (json.uid) !== "undefined") {
			// If job is defined, work unit has been sent
			this.jobCount++;
			console.info("Job " + this.jobCount + " received.");
			this.receiveWorkUnit(json);
			return;
		}

		switch (json.type) {
		case "confirmation":
			console.log("Requesting work unit");
			this.requestWorkUnit();
			break;
		case "retry":
			console.log("Server requests retry.");
			setTimeout(draculaInstance.requestWorkUnit, Math.floor((Math
					.random() * 10000) + 1000));

			break;
		default:
		case "nomore":
			console.log("No work");
			setTimeout(draculaInstance.requestWorkUnit, Math.floor((Math
					.random() * 60000) + 30000));
			break;
		}
	};

	this.isWorkerAvailable = function() {
		return typeof (Worker) !== 'undefined'
	};

	/**
	 * Server Communication
	 */
	this.requestWorkUnit = function() {
		$.getScript(REQUEST_URI + "?r=workunit&callback=" + this.callBack);
	};

	this.receiveWorkUnit = function(json) {
		if (!this.isWorkerAvailable()) {
			var search = new MsSearch(json);
			search.search();
			return;
		}

		this.renfield.postMessage(json);
	};

	this.sendResult = function(resultObject) {
		console.log("Processed in: " + resultObject.processTime);
		var resultString = JSON.stringify(resultObject);

		$.getScript(REQUEST_URI + "?r=result&result=" + resultString
				+ "&callback=" + this.callBack);
	};

	this.sendTerminating = function() {
		this.renfield.terminate(); // terminate the worker...
		$.getScript(REQUEST_URI + "?r=terminate&callback=" + this.callBack);
	};

	this.initialise = function() {
		if (!this.isWorkerAvailable()) {
			$.ajax({
				url : WORKER_URI,
				dataType : 'script',
				async : false
			});

			console.warn("WebWorker not available");
		} else {
			console.info("WebWorker available");
			this.initialiseWorker();
		}
	}

	this.initialiseWorker = function() {
		var BuildWorker = function(foo) {
			var str = foo.toString().match(
					/^\s*function\s*\(\s*\)\s*\{(([\s\S](?!\}$))*[\s\S])/)[1];

			return new Worker(window.URL.createObjectURL(new Blob([ str ], {
				type : 'text/javascript'
			})));
		};

		this.renfield = BuildWorker(function() {
			// will need to go to master
			var SCRIPT_URL = "http://pgb.liv.ac.uk/~andrew/crowdsource-js/src/";
			importScripts(SCRIPT_URL + "MsSearch.js");
		});

		// worker communicates with the main js via JSON strings
		this.renfield.onmessage = function(e) {
			draculaInstance.sendResult(e.data);
		};
	}
}

var draculaClient = new DraculaClient('draculaClient');

$(window).on("beforeunload", function() {
	if (typeof (w) !== "undefined") {
		draculaClient.sendTerminating();
	}
});

draculaClient.initialise();
draculaClient.requestWorkUnit();