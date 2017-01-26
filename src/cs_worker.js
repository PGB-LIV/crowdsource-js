
<<<<<<< HEAD
	
	var MAX_NUMBER_SCORES_RETURNED = 5;			//only retunr the top 5 peptide scores for this workunit
	
	
	var WORKER_INITIALISING = 0;
	var WORKER_AWAITING_WORKUNIT = 1;	
	var WORKER_PROCESSING = 3;
	var WORKER_AWAITING_CONFIRMATION = 4;
	var WORKER_COMPLETE = 5;
=======
//JOHNS DEBUG MODERATOR
	var DEBUG_LEVEL_OFF = 3;
	var DEBUG_LEVEL_HIGH = 2; 
	var DEBUG_LEVEL_MEDIUM = 1;
	var DEBUG_LEVEL_LOW = 0;
	
	var debugMinimumLevel = DEBUG_LEVEL_MEDIUM;	//anything below this priority doesn't get output (set it to HIGH if only critical debug information)
	function Debug(output , priority)	// = DEBUG_LEVEL_HIGH)
	{
		priority = (typeof priority !== 'undefined') ?  priority : DEBUG_LEVEL_HIGH;
		if	(priority >= debugMinimumLevel)
		{					
			console.log(output);
		}
	}
//END

	var WORKER_INITIALISING = 0;
	var WORKER_AWAITING_SPECTRA = 1;	
	var WORKER_AWAITING_PROTEIN = 2;
	var WORKER_PROCESSING = 3;
	var WORKER_AWAITING_CONFIRMATION = 4;
	var WORKER_COMPLETE = 5;
	
	var mySpectrum;
	var myProtein;
	var workerStatus = WORKER_AWAITING_SPECTRA;
	
	

	
	
>>>>>>> c2825958c733843fb6f2abec74df5013487df750
	

	

	
	onmessage = function (event){
		var myObj = JSON.parse(event.data);
		var retString;
		switch (myObj.type)
		{
			case "spectrum":
				retString = '{"type":"acknowledge","what":"spectrum"}';	
				mySpectrum=myObj;
				setStatus(WORKER_AWAITING_PROTEIN);
				break;
			case "protein":
				retString = '{"type":"acknowledge","what":"protein"}';	
				myProtein=myObj;
				setStatus(WORKER_PROCESSING);
				break;
			case "confirmation":
				retString = '{"type":"acknowledge","what":"confirmation"}';	
<<<<<<< HEAD
=======
				setStatus(WORKER_AWAITING_PROTEIN);
>>>>>>> c2825958c733843fb6f2abec74df5013487df750
				break;
			default:
				retString = '{"type":"acknowledge","what":"message"}';	
				break;
		}
		postMessage(retString);
		if (workerStatus === WORKER_PROCESSING)
		{
			doSearch();
		}
	}
	
<<<<<<< HEAD
	//it may be useful to normalise ms2 to 0 - 100 where 100 = maximum intensity.
	function normaliseMs2()
	{
		//warning I actually change mySpectrum!!!
		var _maxIntensity = 0;
		for (var i = 0; i < myWorkUnit.ms2.length; i++)
		{
			if (myWorkUnit.ms2[i].intensity > _maxIntensity)
			{
				_maxIntensity = myWorkUnit.ms2[i].intensity;
			}
		}
		
		for (i = 0; i < myWorkUnit.ms2.length; i++)
		{
			var s = myWorkUnit.ms2[i].intensity;
			s *= 100/_maxIntensity;
			
			myWorkUnit.ms2[i].intensity = s;
		}		
	}

=======
	
>>>>>>> c2825958c733843fb6f2abec74df5013487df750
	function setStatus(newStatus)
	{
		workerStatus = newStatus;
	}
	 
	function doSearch()
	{
		var resultObject = {
			type:"result",
			job:0,
			ms1:0,
			protein: 0,
			score: 0,
			peptides: [],
		}
		
		var pepCandidates =[];	//array of peptide indices (ie their index in protein.peptides[]
		var myB_ions = [];		//array of B ions filled by fragment() upon successful ms1 
		var myY_ions = [];		//array of Y ions filled by fragment() upon successful ms1

		//lets set up the current resultObject properties that we know.
		resultObject.job = myProtein.job;
		resultObject.ms1 = mySpectrum.ms1;
		resultObject.protein = myProtein.protein;
				
		//ms1 search 
		pepCandidates = ms1Search();		//returns with an array of peptide indices (ie their index in protein.peptides[]
		if (pepCandidates.length != 0)		//we have found at least one peptide that matches the pepmass (No PTM as yet) 
		{
			var scores = ms2Search(pepCandidates);			//retuns an array of objects [{cumIntensity:0,b_ions:0,y_ions:0}];
			for (i = 0; i < pepCandidates.length; i++)
			{
				var peptObj = {id:0,score:0,b_count:0,y_count:0};
				peptObj.id=myProtein.peptides[pepCandidates[i]].id;
				peptObj.score = scores[i].cumIntensity;
				peptObj.b_count = scores[i].b_ions;
				peptObj.y_count=scores[i].y_ions;
			
				resultObject.peptides.push(peptObj);
			}
		}
		var retString = JSON.stringify(resultObject);
		postMessage(retString);
	}		
	
	function ms1Search(){
		var candidates = [];
		var myPepmass = mySpectrum.pepmass*mySpectrum.charge - (mySpectrum.charge * 1.007276);		//adjusted pepmass = ms1 pepmass * charge - charge*H+; 
		var minPepmass = myPepmass - (myPepmass*.000010);
		var maxPepmass = myPepmass + (myPepmass*.000010);
		//Debug("number of peptides checked against = " +myProtein.peptides.length);
		for (i = 0; i < myProtein.peptides.length; i++)
		{
			var m =myProtein.peptides[i].mass;
			if (m >= minPepmass && m <= maxPepmass)
			{
				candidates.push(i);
			}
		}
		return candidates;		//simple array of indexes of candidate peptides 
	}	
	
	
	
	function ms2Search(candidates)
	{
		
<<<<<<< HEAD
				
			//lets set up the current default resultObject properties that we know.
		resultObject.job = myWorkUnit.job;
		resultObject.workunit = myWorkUnit.id;
		for (p = 0; p < myWorkUnit.peptides.length; p++)
			{
				var scoreObj = {id:0, score:0};
				bscore = 0;
				yscore = 0;
				bcount = 0;
				ycount = 0;
				
				fragment(myWorkUnit.peptides[p].structure);			//flush out myB_ions & myY_ions array.
				
				for (b = 0; b < myB_ions.length; b++)
=======
		var scores = [];
		Debug("number of candidates = "+candidates.length);
		for (j = 0; j < candidates.length; j++)
		{
			var scoreObj = {cumIntensity:0,b_ions:0,y_ions:0};
			score = 0;
			bcount = 0;
			ycount = 0;
			
			Debug("peptide of protein = "+candidates[j]);
			fragment(candidates[j]);		//flush out myB_ions & myY_ions array.
											//score+ = intensity of a match (+-0.2Daltons)
											//check B ions against ms2 records;
			Debug("B ions = "+myB_ions,DEBUG_LEVEL_LOW);
			Debug("Y ions = "+myY_ions.DEBUG_LEVEL_LOW);
						
			for (b = 0; b < myB_ions.length; b++)
			{
				for (m = 0; m < mySpectrum.ms2.length; m++)
>>>>>>> c2825958c733843fb6f2abec74df5013487df750
				{
					if (Math.abs(myB_ions[b]-mySpectrum.ms2[m]['mz']) <= 0.2)
					{
<<<<<<< HEAD
						if (Math.abs(myB_ions[b]-myWorkUnit.ms2[m]['mz']) <= 0.2)
						{
							bcount++;
							bscore += myWorkUnit.ms2[m]['intensity'];		
						}
=======
						Debug("B ion match");
						bcount++;
						score += mySpectrum.ms2[m]['intensity'];		
>>>>>>> c2825958c733843fb6f2abec74df5013487df750
					}
				}
			}
											//check Y ions against ms2 records;								
			for (y = 0; y < myY_ions.length; y++)
			{
				for (m = 0; m < mySpectrum.ms2.length; m++)
				{
					if (Math.abs(myY_ions[y]-mySpectrum.ms2[m]['mz']) <= 0.2)
					{
<<<<<<< HEAD
						if (Math.abs(myY_ions[y]-myWorkUnit.ms2[m]['mz']) <= 0.2)
						{
							ycount++;
							yscore += myWorkUnit.ms2[m]['intensity'];		
						}
					}
				}
				
				//important bit.... 				//score = (bscore+yscore)*(bcount+ycount);		//we will need to investigate the best scoring formula 
				scoreObj.score = (bscore+yscore)*(bcount)*(ycount);
				scoreObj.id = myWorkUnit.peptides[p].id;
				allPeptideScores.push(scoreObj);
			}
	
		//Only wish to return the MAX_NUMBER_SCORES_RETURNED	
		//Order array by score 
			
			allPeptideScores.sort(function(a,b){return b.score - a.score});
			var numRes = MAX_NUMBER_SCORES_RETURNED;
			if (myWorkUnit.peptides.length < numRes)
			{
				numRes = myWorkUnit.peptides.length;
			}
			
			for (r = 0; r < numRes; r++)
			{
				resultObject.peptides.push(allPeptideScores[r])
			}
			var retString = JSON.stringify(resultObject);
			console.log("job = "+resultObject.workunit);
			postMessage(retString);
=======
						Debug("Y ion match");
						ycount++;
						score += mySpectrum.ms2[m]['intensity'];		
					}
				}
			}
			Debug("score ="+score);
			scoreObj.cumIntensity = score;
			scoreObj.b_ions=bcount;
			scoreObj.y_ions=ycount;
			scores.push(scoreObj);
		}
		return scores;
>>>>>>> c2825958c733843fb6f2abec74df5013487df750
	}
		
	

			
	var AAMass = {A:71.037114,
				R:156.101111,
				N:114.042927,
				D:115.026943,
				C:103.009185,
				E:129.042593,
				Q:128.058578,
				G:57.021464,
				H:137.058912,
				I:113.084064,
				L:113.084064,
				K:128.094963,
				M:131.040485,
				F:147.068414,
				P:97.052764,
				S:87.032028,
				T:101.047679,
				U:150.95363,
				W:186.079313,
				Y:163.063329,
				V:99.068414
				};
				
	
	//creates arrays of B and Y ions in myB and myY of peptide[index] held in myProtein;
	function fragment(index)				 
	{
<<<<<<< HEAD
		var myPeptide = struc;		//myProtein.peptides[index]['structure'];
=======
		
		var myPeptide = myProtein.peptides[index]['structure'];
		Debug("Fragmenting "+myPeptide+ " and id = "+myProtein.peptides[index]["id"]);
>>>>>>> c2825958c733843fb6f2abec74df5013487df750
		var cm = 0;
		myB_ions = [];
		for (i = 0; i < myPeptide.length; i++)
		{
			var aa = myPeptide.charAt(i);	
			var m = AAMass[aa]; 
			if (i === 0){	//first time through? -0H + water
				m += 1.007276;		//add a hydrogen
			}
			cm+=m;
			myB_ions[i]= cm.toFixed(6);
		}
		cm = 0;
		myY_ions=[];
		for (var i = myPeptide.length-1; i >= 0; i--)
		{
<<<<<<< HEAD
			var aay = myPeptide.charAt(i);
			modmass = checkforPTM(aay,'y');
			var my = AAMass[aay]; 
			if (i === (myPeptide.length-1)){	//first time through? H + water
				my +=  18.010565+1.007276;		
			}
			my+=modmass;
			cm+=my;
=======
			var aa = myPeptide.charAt(i);
			var m = AAMass[aa]; 
			if (i == (myPeptide.length-1)){	//first time through? H + water
				m +=  18.010565+1.007276;		
			}
			cm+=m;
>>>>>>> c2825958c733843fb6f2abec74df5013487df750
			myY_ions[i]= cm.toFixed(6);
		}
	}
	

<<<<<<< HEAD
	function checkforPTM(aap,ion)
	{
		var masstoadd = 0;
		for (var m = 0; m < myWorkUnit.mods.length; m++)
		{
			if (myWorkUnit.mods[m]['loc']==aap)
			{
				masstoadd+=myWorkUnit.mods[m]['modmass'];
			}
		}
		return masstoadd;
	}
=======

>>>>>>> c2825958c733843fb6f2abec74df5013487df750
		
		
		
		
		
		
		
		
	
		
		