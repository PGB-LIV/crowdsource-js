

	var PROVIDER_PHP = "http://pgb.liv.ac.uk/~johnheap/crowdsource-server/src/public_html/testscheduler.php";		//the php that handles data requests and accepts results (will need to go to master)

	var myWorker;			//WebWorker instance		
	var myWorkUnit;
	
	var BuildWorker = function(foo){
	var str = foo.toString()
             .match(/^\s*function\s*\(\s*\)\s*\{(([\s\S](?!\}$))*[\s\S])/)[1];
	return  new Worker(window.URL.createObjectURL(
                      new Blob([str],{type:'text/javascript'})));
	}

	
	initialiseWorker();
	requestWorkUnit();
			
	function initialiseWorker()
	{
		myWorker = BuildWorker(function()
		{	
			var WORKER_NAME = "http://pgb.liv.ac.uk/~johnheap/crowdsource-server/src/public_html/javascript/cs_worker.js";			//name of worker will need to goto master
			importScripts(WORKER_NAME);
		});
	}
		
	function requestWorkUnit()
	{
		$.getScript(PROVIDER_PHP+"?r=workunit");
	}
	
	function receiveWorkUnit(json)
	{
		myWorkUnit = JSON.stringify(json);
		myWorker.postMessage(myWorkUnit);		//on receiving workUnit it will get to work	
	}
	
	
	function sendResult(resultString)
	{
		$.getScript(PROVIDER_PHP+"?r=result&result="+resultString);
	}
	
	function sendTerminating()
	{
		$.getScript(PROVIDER_PHP+"?r=terminate");
	}
	
	
	//this function is the P of the JSONP response from server eg "parseResult(Object)"
	function parseResult(json)
	{
		switch (json.type)
		{
			case "workunit":

				receiveWorkUnit(json);
				break;
				
			case "confirmation":
				requestWorkUnit();
				break;
			
			case "message":
				break;
			case "nomore":
				break;
			default:
				break;
		}
		json=null;
	}

	myWorker.onmessage = function(e)		//worker communicates with the main js via JSON strings
	{
		var workerResponse = JSON.parse(e.data);
		
		switch (workerResponse.type)
		{
			case "acknowledge":
			switch (workerResponse.what)
			{
				case "workunit": 
					break;
				case "confirmation":		
					break;
			}
			break;
			
			case "result":
				sendResult(JSON.stringify(workerResponse));
				break;
			default:
				break;
		}
	}

	//terminate session event. Also occurs on Refresh.
	$(window).on("beforeunload", function() {
		if (typeof(w) !== "undefined")
		{
			myWorker.terminate();						//terminate the worker...
			sendTerminating();
		}
	});			
	

				
			
				
	
	


