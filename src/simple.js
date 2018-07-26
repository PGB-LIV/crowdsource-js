

var TEST_PHP = "http://pgb.liv.ac.uk/~johnheap/crowdsource-server/src/public_html/memorytest.php";


//var t = setInterval(trySimple,200);
	trySimple();
	
	
function trySimple()
{
//	console.log("inserting script");
	$.getScript(TEST_PHP);
}


function parseResult(json)
{
	setTimeout(trySimple, 200);	
}


