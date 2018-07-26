"use strict";
/**
 * Copyright 2016 University of Liverpool
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
var PHOSPHORYLATION = 21;				//global equate. (may change dependent on Uniprot!!!!)
var PHOSPHO_LOSS = 79.966;  			//Phosphorylation is a special case hence special globals ....
										//in actual fact probably all modifications require special case scoring due to different proton affinities etc.
										//but phophorylation definitely does so rather than have a 500 entry jump table /switch statement we'll try this first.

	//Redoing strategy as this isn't really phase 2 but phase1 & 2
	//Most peptide candidates may not be present... 
	//Candidates chosen precursor == peptide mass + (n *PHOSPHO) where n = number of modifications 
	
	function doSecondPhaseSearch(myWorkUnit)
	{
		var BAIL_BELOW_SCORE = 1;		//do not continue peptide evaluation if scores less than this
		var allPeptideScores = [];		//made up of {peptide id, ionsmatched, score  mods:[{ id,position[] }]}
		
		for (var p = 0; p < myWorkUnit.peptides.length; p++){
		//for (var p = 0; p < 2; p++){
			var scoreObj = {};
			var currPeptide = myWorkUnit.peptides[p];
			var resObj = {id:currPeptide.id,ionsMatched:0,score:0, mods:[{}]};
			
			for (var mod = 0; mod < currPeptide.mods.length; mod++){	//we may have more than one possible modification per peptide as some point!!!
				var currMod = currPeptide.mods[mod];
				//FUDGE
				if (currMod.id === 21){
					//currMod.residues = "STY";
				}
				//END FUDGE
				var possLocs  = createIndexofPossibleLocations(currPeptide.sequence,currMod.residues);		//possLocs an array of possible modification locations+1
				possLocs.sort(function(a, b){return a - b});			//sort the array 0 -> seqlen
				//console.log("possLocs:"+JSON.stringify(possLocs));
				var confirmedPos = [];									//holds the confirmed locations of a mod (via matching Ions) 
				var currScoreObj = {modPos:[],ionsMatched:0,score:0}; 	//holds the best so far score for this peptide modification )positions held in modpos[]
																		
				for (var num = 0; num <= currMod.num; num++){											//as some mods fall off we are going to check for 0 mods to num mods	
					var subIonsets = getIonsetForMod(currPeptide,currMod,possLocs,confirmedPos,num);	//an array of ionsets to score containing num number of mods
					for (var s = 0; s < subIonsets.length;s++){	
						subIonsets[s]= matchSpectraWithIonSet(myWorkUnit.fragments, subIonsets[s]);		//for each ionset with log matches with ms2 fragments
						scoreObj = scoreIonset(subIonsets[s]);											//score the ionset (matched ions and ion ladders)
						confirmedPos = confirmMods(subIonsets[s],confirmedPos);							//see if we have any matches that confirm a mod. (bit unsophisticated could do better as modmass continues along ion range)
						
						
						if (scoreObj.score > currScoreObj.score){
							currScoreObj.modPos = subIonsets[s].modPos.slice();
							//for (var i = 0; i < subIonsets[s].modPos.length;i++){
							//	currScoreObj.modPos[i] = subIonsets[s].modPos[i];
							//}
						}
						currScoreObj.score =(scoreObj.score > currScoreObj.score)?scoreObj.score:currScoreObj.score;  //take the best score so far
						currScoreObj.ionsMatched =(scoreObj.ionsMatched > currScoreObj.ionsMatched)?scoreObj.ionsMatched:currScoreObj.ionsMatched;  //take the best ionsmatched so far
					}
					if (currScoreObj.score < BAIL_BELOW_SCORE || subIonsets.length ===0){
						//console.log ("Doing no more with Current peptide id:"+currPeptide.id);
						break;		//bail if score to small or ionsets.len = 0					
					}	
				}
				
				resObj.score += currScoreObj.score;  
				resObj.ionsMatched = (currScoreObj.ionsMatched>resObj.ionsMatched)?currScoreObj.ionsMatched:resObj.ionsMatched;
				var resMod = {id:currMod.id, position:currScoreObj.modPos};	
				resMod.position = uniqueArrayConcat(resMod.position, confirmedPos, resMod.position.length-currMod.num);		
				
				if (confirmedPos.length > 0){
					console.log("Peptide:"+currPeptide.id+", "+currPeptide.sequence+"\n" +
							JSON.stringify(resObj)+"\n---------");
//					console.log("Peptide:"+currPeptide.id+", "+currPeptide.sequence+", "+JSON.stringify(possLocs)+", ["+currMod.num+"]\n" +
//							"Score:"+resObj.score+" Match:"+resObj.ionsMatched+" loc:"+JSON.stringify(resMod.position)+"\n---------");
				}
				
				/*
				if 	(resMod.position.length < currMod.num){			//need more mod positions....
					//if score is reasonable... if num = number of plocs... then we could infer their positions ..  but not increase score
					if (resObj.score > 3){
						possLocs  = createIndexofPossibleLocations(currPeptide.sequence,currMod.residues);		//get it again becasue I 've screwed it.
						if (possLocs.length === currMod.num){
							resMod.position = possLocs.slice();
						}
					}
				}
				*/
					
								
				//For the moment I have to bulk out the position field to be equal to number of mods (if >=200 then I haven't found/cannot confirm a position)
				if 	(resMod.position.length < currMod.num){								//(confirmedPos.length < currMod.num){
					var n = currMod.num - resMod.position.length; //confirmedPos.length;
					for (var t =0; t < n; t++){
						resMod.position.push(200+t);
					}
				}
				resObj.mods[mod]=resMod;
			//	console.log("resObj:"+JSON.stringify(resObj));
			}
			allPeptideScores.push(resObj);
		}
		
		var num2Res = 50;		//maximum number of results from phase2
		var resultObject= {job:myWorkUnit.job,precursor:myWorkUnit.precursor,peptides:[]};		//what is sent back as retString;

		allPeptideScores.sort(function(a,b){return b.score - a.score});		
		
		if (myWorkUnit.peptides.length < num2Res) {
			num2Res = myWorkUnit.peptides.length;
		}
		for (var r = 0; r < num2Res; r++) {
			resultObject.peptides.push(allPeptideScores[r])
		}
		var retString = JSON.stringify(resultObject);
		console.log("The return string = "+retString);
		console.log("pre = "+resultObject.precursor);
		if (resultObject.precursor < 1001){			//for testing
		//	postMessage(retString);
		}
		postMessage(retString);
		return resultObject; //for test
	}
	
	
	
	
	
	function getIonsetForMod(myPeptide,myMod,poslocs,conpos,num)
	{	
		var ionSetArray = [];		//array ionsets to return
		var combArray = [];			//combination array when dealing with more than one number of modifications of one type 
		var modMass = myMod.mass	//shorthand;
		var MUCH_TOO_MUCH = 3500;	//Arbitrary value at which point I can't be arsed building so many ionsets.
		
		
		//search for none modified peptide - may well be here (Peaks PTM) and gives a basal score for peptide.
		if (num === 0){
			
			var ions = getIonsFromArray(myPeptide,[],0);			//mypeptide,locs[],mass
			ionSetArray.push(ions);									//single ionset returned
			return ionSetArray;
		}
		//search for occurrence of a single modification at any of the possible possible locations
		//the others may have fallen off. Finding a good ion match can reduce size of future search
		if (num === 1){
			for (var pl = 0; pl < poslocs.length; pl++){
				var ions = getIonsFromArray(myPeptide,[poslocs[pl]],modMass);		//an array of one!
				ionSetArray.push(ions);
			}
			return ionSetArray;
		}
	
		//we must have got a hit (or two) as we've continued
	
		if (conpos.length >= num){			//first check that we haven't found all we can. (as if!)
			return ionSetArray;
		}
		//if we have any confirmed positions, take them out of poslocs as we don't need them for the combinations.
		for (var cp = 0; cp < conpos.length; cp++){
			var p = conpos[cp];
			for (var pl = 0; pl < poslocs.length; pl++){
				if (poslocs[pl] === p){
					poslocs.splice(pl,1);	//remove it.	
				}
			}		
		}
		
		//do a calculation first to see if do-able combinations rapidly get out of hand.
		//number of combinations = n!/(k!(n-k)!)
		if (getCombinations(poslocs.length,num) > MUCH_TOO_MUCH){
			console.log("Bailing Array too large for me. Combinations of "+poslocs.length+","+num+"="+getCombinations(poslocs.length,num));
			return ionSetArray;
		}
		
		combArray = createArrayPossibleCombinations(poslocs, poslocs.length,num-conpos.length);
		
		//need to add confirmed positions (conpos) back into the mix. 
		for (var c = 0; c < combArray.length; c++){	
			var mlocs = conpos.concat(combArray[c]);		//mlocs = confirmed mods + combination to try
			var ions = getIonsFromArray(myPeptide,mlocs,modMass);
			ionSetArray.push(ions);
		}
		return ionSetArray;
	}
	
	function getIonsFromArray(myPeptide,plocs,vmodmass)		//looking for more than one modification at a time
	{
		var ionSet = new Ionset();			//{modPos:[],bIons:[],yIons:[]};		//bIons = [{mass, match, intensity}]		
		var acid = '';						//amino acid character
		var sequence = myPeptide.sequence; 	//'ABCEF'
		
		ionSet.modPos = plocs.slice();		//initialise modpos;


		var cumulativeMass = 1.007276;				//cumulative mass for b ions start with Water - OH which is H
		cumulativeMass+=checkforFixedPTM('[');		//is a fixed mod an amine terminus?
		for (var b = 0; b < (sequence.length); b++)
		{
			var ionObj = new Ion();		//{mass:0,match:0,intensity:0,deltaM:0,modFlag:false};
			acid = sequence.charAt(b);
			cumulativeMass+= g_AAmass[acid]+checkforFixedPTM(acid);
			for (var pl = 0; pl < plocs.length; pl++){
				if (b === (plocs[pl]-1)){
					cumulativeMass+=vmodmass;
					ionObj.modFlag = true;
					
				}
			}
			ionObj.mass = cumulativeMass;
			
			ionSet.bIons.push(ionObj);
		}
		
		cumulativeMass = 18.010565+1.007276;			//start with H + water
		cumulativeMass+= checkforFixedPTM(']');			//is a fixed mod a carboxyl end?
		for (var y = sequence.length-1; y >= 0; y--)
		{
			var ionObj = new Ion();	//{mass:0,match:0,intensity:0,deltaM:0,modFlag:false};
			acid = sequence.charAt(y);
			cumulativeMass+=g_AAmass[acid]+checkforFixedPTM(acid);
			for (var pl = 0; pl < plocs.length; pl++){
				if (y === (plocs[pl]-1)){
					cumulativeMass+=vmodmass;
					ionObj.modFlag = true;
					
				}
			}
			ionObj.mass = cumulativeMass;
			
			ionSet.yIons.push(ionObj);
		}
		return ionSet;
	}
	
	
	//modifies an array of confirmed modification locations (conpos) should we find a confirmation (both modFlag and match) in this ionset.  
	function confirmMods(ionset,conpos)
	{
		var contemp = [];
		var modPos = -1;
		for (var b = 0; b < ionset.bIons.length; b++){
			if (ionset.bIons[b].modFlag){
				modPos = b+1;
			}
			if (modPos !== -1 && ionset.bIons[b].match){		//all matches after a modification ion result in a modpos being stored in contemp[]
				contemp.push(modPos);
			}
		}
		modPos = -1;
		for (var y = 0; y < ionset.yIons.length; y++){
			if (ionset.yIons[y].modFlag){
				modPos=ionset.yIons.length-1-y+1;
			}
			if (modPos !== -1 && ionset.yIons[y].match){
				contemp.push(modPos);
			}
		}
		for (var ct = 0; ct < contemp.length; ct++){		//check new locations held in contemp[] to see if already in conpos 
			var p = contemp[ct];
			var pmatch = false;
			for (var cp = 0; cp < conpos.length; cp++){
				if (conpos[cp] === p){
					pmatch =true;							//already in conpos so no need to add this 
				}
			}
			if (!pmatch){
				conpos.push(p)								//new one found in this ions set so add to end;
			}
		}
		return conpos;
	}
	/*function confirmMods(ionset,conpos)
	{
		var contemp = [];
		for (var b = 0; b < ionset.bIons.length; b++){
			if (ionset.bIons[b].modFlag && ionset.bIons[b].match)
			{
				contemp.push(b+1);		//found a confirmed modification location
			}
		}
		for (var y = 0; y < ionset.yIons.length; y++){
			if (ionset.yIons[y].modFlag && ionset.yIons[y].match){
				contemp.push (ionset.yIons.length-1-y+1);	//found a confirmed modificaton location
			}
		}
		for (var ct = 0; ct < contemp.length; ct++){		//check new locations held in contemp[] to see if already in conpos 
			var p = contemp[ct];
			var pmatch = false;
			for (var cp = 0; cp < conpos.length; cp++){
				if (conpos[cp] === p){
					pmatch =true;							//already in conpos so no need to add this 
				}
			}
			if (!pmatch){
				conpos.push(p)								//new one found in this ions set so add to end;
			}
		}
		return conpos;
	}*/
	
	
	//returns integer array of all possible locations of mod +1 
	//NB position = sequence position + 1 ... 0 = '[' 1 to len = = to len-1 and length+1 = ] ] 
	//input sequence (eg "PEPTIDE") and residues (eg "CKV");
	function createIndexofPossibleLocations(sequence, residues){
		var posloc = [];		//array of possible locations for this mod
		var pos=0;
		var res='';
		//special case for residues '[' annd ']'
		for (var r = 0; r < residues.length; r++){
			res = residues[r];
			if (res === '['){
				pos = 0;
				posloc.push[pos];
			}else if (res === ']'){
				pos = residues.length+1;
				posloc.push[pos];
			}else{
				for (pos = 0; pos < sequence.length; pos++){
					if (sequence[pos] === res){
						posloc.push(pos+1);				//using modified location 1 = position zero
					}
				}
			}
		}
		return posloc;
	}
	
	
	
	
	//--------------------------------------------------------	
	//adds unique parts of second to first until maxlen is reached
	function uniqueArrayConcat(first,second,maxlen)
	{
		if (first.length >= maxlen){
			return first;				
		}
		
		for (var s = 0; s < second.length; s++){		//check new locations held in contemp[] to see if already in conpos 
			var val = second[s];
			var valmatch = false;
			for (var f = 0; f < first.length; f++){
				if (first[f] === val){
					valmatch =true;							//already in conpos so no need to add this 
				}
			}
			if (!valmatch){
				if (first.length < maxlen){
					first.push(val);								//new one found in this ions set so add to end;
				}
			}
		}
		return first;
	}
	
	//-----------------------------------------------------------------------------
	// Combination & Factorial code
	
	function createArrayPossibleCombinations(locArray, n, k) {
		var result = [];
		var values = [];
		for (var i = 0; i <= n; i++) {
		  values[i] = i;
		}
			// initialize permutations
		var perm = [];
		for (var i = 0; i < n; i++) {
		  if (i < k) {
		    perm[i] = 1;
		  } else {
		    perm[i] = 0;
		  }
		}
		perm.sort();

		whileloop:
		  while (true) {
		    // save subresult
		    var subresult = [];
		    for (var i = 0; i < n; i++) {
		      if (perm[i] === 1) {
		        subresult.push(values[i]);
		      }
		    }
		    result.push(subresult);
		    // get next permutation
		    for (var i = n - 1; i > 0; i--) {
		      if (perm[i - 1] === 1) {
		        continue;
		      }
		      if (perm[i] === 1) {
		        perm[i - 1] = 1;
		        perm[i] = 0;
		        perm = perm.slice(0, i).concat(perm.slice(i).sort())
		        continue whileloop;
		      }
		    }
		    // no additional permutations exist
		    break whileloop;
		  }
		  //result = array of integers that correspond to index in location array 
		  for (i = 0; i < result.length; i++){
			  for (var j = 0; j < k; j++){
				  result[i][j] = locArray[result[i][j]];  //translate to actual locations.
			  }
		  }
	  return result;			
	}
	function getCombinations(n,k)
	{
		var nfac= getFactorial(n);
		var kfac = getFactorial(k);
		var nsubkfac = getFactorial(n-k);
		return nfac/(kfac*nsubkfac);
	}
	function getFactorial(f)
	{
		var r = 1
		for (;f >0;r*=f,f--){}
		return r;
	}
		
	
	
		

	