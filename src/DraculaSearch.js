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

	this.isNoWorkerAllowed = false;
	var REQUEST_URI = 'http://pgb.liv.ac.uk:1260/work';
	var WORKER_URI = 'http://pgb.liv.ac.uk:1260/script/MsSearch.js';
	var VERSION = 'INSERT_BUILD_VERSION';
	var draculaInstance = this;

	this.callBack = callBackInstance + '.parseResult';

	this.renfield = null;

	this.jobCount = 0;

	/**
	 * JSONP method
	 */
	this.parseResult = function(json) {
		if (typeof (json.uid) !== 'undefined') {
			// If job is defined, work unit has been sent
			this.jobCount++;
			console.info('Job ' + this.jobCount + ' received.');
			this.receiveWorkUnit(json);
			return;
		}

		switch (json.type) {
		case 'confirmation':
			console.log('Requesting work unit');
			this.requestWorkUnit();
			break;
		case 'retry':
			console.log('Server requests retry.');
			setTimeout(draculaInstance.requestWorkUnit, Math.floor((Math
					.random() * 10000) + 1000));

			break;
		case 'nomore':
			console.log('No work');
			setTimeout(draculaInstance.requestWorkUnit, Math.floor((Math
					.random() * 60000) + 30000));
			break;
		case 'BadVersion':
		default:
			throw 'Client outdated.';
			break;
		}
	};

	this.isWorkerAvailable = function() {
		return typeof (Worker) !== 'undefined';
	};

	/**
	 * Server Communication
	 */
	this.requestWorkUnit = function() {
		draculaInstance.getScript(REQUEST_URI + '?v=' + VERSION
				+ '&r=workunit&callback=' + draculaInstance.callBack);
	};

	this.receiveWorkUnit = function(json) {
		if (!this.isWorkerAvailable() && this.isNoWorkerAllowed) {
			var search = new MsSearch(json);
			search.search();
			return;
		}

		this.renfield.postMessage(json);
	};

	this.sendResult = function(resultObject) {
		console.log('Processed in: ' + resultObject.processTime);
		var resultString = JSON.stringify(resultObject);

		this.getScript(REQUEST_URI + '?v=' + VERSION + '&r=result&result='
				+ resultString + '&callback=' + this.callBack);
	};

	this.sendTerminating = function() {
		this.renfield.terminate(); // terminate the worker...
		this.getScript(REQUEST_URI + '?r=terminate&callback=' + this.callBack);
	};

	this.initialise = function() {
		if (!this.isWorkerAvailable()) {
			this.getScript(WORKER_URI);

			console.warn('WebWorker not available');
		} else {
			console.info('WebWorker available');
			this.initialiseWorker();
		}
	};

	this.initialiseWorker = function() {
		var BuildWorker = function(importFunc) {
			var strImportFunc = 'var importScript = \'' + WORKER_URI + '\';';
			strImportFunc += importFunc.toString().match(
					/^\s*function\s*\(\s*\)\s*\{(([\s\S](?!\}$))*[\s\S])/)[1];

			return new Worker(window.URL.createObjectURL(new Blob(
					[ strImportFunc ], {
						type : 'text/javascript'
					})));
		};

		this.renfield = BuildWorker(function() {
			importScripts(importScript);
		});

		this.renfield.onmessage = function(e) {
			draculaInstance.sendResult(e.data);
		};
	};

	// https://stackoverflow.com/a/28002292/702192
	this.getScript = function(source, callback) {
		var script = document.createElement('script');
		var prior = document.getElementsByTagName('script')[0];
		script.async = 1;

		script.onload = script.onreadystatechange = function(_, isAbort) {
			if (isAbort || !script.readyState
					|| /loaded|complete/.test(script.readyState)) {
				script.onload = script.onreadystatechange = null;
				script = null;

				if (!isAbort && callback) {
					setTimeout(callback, 0);
				}
			}
		};
		script.onerror = function() {
			console.warn("Connection failed. Retrying.");
			setTimeout(draculaInstance.requestWorkUnit, Math.floor((Math
					.random() * 10000) + 1000));
		};

		script.src = source;
		prior.parentNode.insertBefore(script, prior);
	};
}

var draculaClient = new DraculaClient('draculaClient');

draculaClient.initialise();
draculaClient.requestWorkUnit();
