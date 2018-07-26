"use strict";
var WORKER_LOCATION = "http://pgb.liv.ac.uk/~johnheap/crowdsource-server/src/public_html/javascript/";  //will need to go to master
var WORKER1_NAME = "cs_worker.js";			//name of worker javascript
var WORKER2_NAME = "firstphase_worker.js";
var WORKER3_NAME = "secondphase_worker.js";

importScripts(WORKER_LOCATION+WORKER1_NAME);
importScripts(WORKER_LOCATION+WORKER2_NAME);
importScripts(WORKER_LOCATION+"thirdphase_worker.js");

