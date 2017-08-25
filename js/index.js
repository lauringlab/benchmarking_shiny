var config = {
    Titer : { value: 5, min: 3, max: 5, step: 1 },
    Phred: { value: 30, min: 30, max: 42, step: 1 },
    MapQ: { value: 0, min: 0, max: 42, step: 1 },
    ReadPosMin:{value:0, min:0, max :62, step:1},
    ReadPosMax:{value:125,min:62,max:125,step:1}
};


$('.slider').each(function(){
    var thisConfig = config[this.id],
        $thisVal = $('#'+this.id+'_val');
    thisConfig.slide = function(event, ui) {
        $thisVal.html(ui.value);
        go();
    };
    $(this).slider(thisConfig);
    $thisVal.html(thisConfig.value);
});

function getVals() {
    var vals = {};
    $('.slider').each(function(){
        vals[this.id] = +($('#'+this.id+'_val').text());
    });
    return vals;
}




function FileHelper()
{}
{
    FileHelper.readStringFromFileAtPath = function(pathOfFileToReadFrom)
    {
        var request = new XMLHttpRequest();
        request.open("GET", pathOfFileToReadFrom, false);
        request.send(null);
        var returnValue = request.responseText;

        return returnValue;
    }
}
function size(ar){
    var row_count = ar.length;
    var row_sizes = []
    for(var i=0;i<row_count;i++){
        row_sizes.push(ar[i].length)
    }
    return [row_count, Math.min.apply(null, row_sizes)]
}

function getdata(file){
	var csv = FileHelper.readStringFromFileAtPath ( file ),
	 data = $.csv.toArrays(csv),
	 header = data[0],
	 allMapQ = [],
	 allGc = [],
	 allPhred = [],
	 allCategory =[],
   allPos = [];
   allPval=[];
   allExp = [];
	for(i=0;i<header.length;i++){
		if(header[i]==="MapQ"){
		var mapqIndex = i;
		}else if( header[i]==="Phred"){
		var phredIndex=i;
		}else if (header[i]==="category"){
		var categoryIndex=i;
		}else if (header[i]==="gc"){
		var gcIndex=i;
		}else if (header[i]==="p.val"){
		var pvalIndex=i;
  } else if (header[i]=="Read_pos"){
    var readIndex=i;
  }else if (header[i]=="exp.freq"){
    var freqIndex=i;
  }
	};
  //console.log(readIndex);
	//indexes=[mapqIndex,	phredIndex,categoryIndex,gcIndex,pvalIndex];
	//console.log(indexes);
	dim = size(data)[0];
	for(i=1;i<dim;i++){ // Make sure this includes the last one
		if( data[i][pvalIndex]<0.01){
		allMapQ.push(data[i][mapqIndex]);
		allGc.push(data[i][gcIndex]);
		allPhred.push(data[i][phredIndex]);
		allCategory.push(data[i][categoryIndex]);
    allPos.push(data[i][readIndex]);
    allPval.push(data[i][pvalIndex]);
    allExp.push(data[i][freqIndex]);
		}
	}
	relevantData=[allGc,allCategory,allMapQ,allPhred,allPos,allPval,allExp];
	//console.log(relevantData[0])
	return(relevantData);
};
var plotData = getdata("https://media.githubusercontent.com/media/lauringlab/benchmarking_shiny/master/20_mut_data_set/processed_data/no.dups.bonferroni.two.sided.sum.csv");


function selectData(titer,mapqCut,phredCut,readCutMin,readCutMax,data) {
	var gc= data[0],
	category = data[1],
	mapq = data[2],
	phred = data[3],
  pos = data[4],
  pval = data[5],
  exp = data[6],
	currentTpMapQ = [],
	currentTpPhred = [],
	currentFpMapQ=[],
	currentFpPhred=[],
  currentFpPos = [],
  currentTpPos = [],
  currentTpPval = [],
  currentFpPval=[],
  currentTpExp=[],
  currentFpExp = [];
  for(i=0;i<gc.length;i++){


		if ( parseInt(gc[i]) === titer && category[i]==="TRUE" ) {

			currentTpMapQ.push(mapq[i]);
			currentTpPhred.push(phred[i]);
      if(mapq[i]>=mapqCut && phred[i]>= phredCut){
        currentTpPos.push(pos[i]);
        if(pos[i]<=readCutMax && pos[i]>=readCutMin){
          currentTpPval.push(pval[i]);
          currentTpExp.push(exp[i]);
        }
      }

		}else if (parseInt(gc[i]) === titer && category[i]==="FALSE" ) {
			currentFpMapQ.push(mapq[i]);
			currentFpPhred.push(phred[i]);
      if(mapq[i]>=mapqCut && phred[i]>= phredCut){
        currentFpPos.push(pos[i]);
        if(pos[i]<=readCutMax && pos[i]>=readCutMin){
          currentFpPval.push(pval[i]);
          currentFpExp.push(exp[i]);
        }
      }
		}
    }
    //console.log(currentTpPval.length);
    //console.log(currentFpPval.length);


  return  [currentFpPhred, currentFpMapQ ,currentTpPhred,currentTpMapQ,currentFpPos,currentTpPos,currentFpPval,currentTpPval,currentFpExp,currentTpExp];

}




function bin(x,bin){
  return Math.round(x/bin)*bin;
}
function binReads(data,binSize){
  var binedData = [];
  for(i=0;i<data.length;i++){
    binedData.push(bin(data[i],binSize));
  }
  return binedData;
}
Array.prototype.unique = function()
{
	var n = {},r=[];
	for(var i = 0; i < this.length; i++)
	{
		if (!n[this[i]])
		{
			n[this[i]] = true;
			r.push(this[i]);
		}
	}
	return r;
}

function count(arr, vals) { //vals is an array of the unique entries in arr
    totalCounts=[];
    for(i=0;i<vals.length;i++){
      var counter=0;
      for(j=0;j<arr.length;j++){
        if(arr[j]===vals[i]){

          counter++;
          //arr.splice(j,1)
        }
        //console.log("moving on")
      }
      totalCounts.push(counter);
    }
    return totalCounts
}

function makeSteps(x,y){
  var newX = [x[0]], //need the first points
      newY=[y[0]];
  for(i=1;i<x.length;i++){
    newX.push(x[i]);
    newX.push(x[i]);   //adds another xi at yi-1
    newY.push(y[i-1]);
    newY.push(y[i]);

  }
return [newX,newY]; //
}

function calcROC(tpPval,fpPval){
 var tpPvalSort = tpPval.sort(function(a, b){return a-b}),
    fdr=[0],//everything starts at 0,0
    possible_fp = 40104,
    possible_tp = 20,
    tdr = [0];
//console.log(tpPvalSort);
  for(i=0; i<tpPvalSort.length;i++){
    var tpc = 0,
        fpc =0;
        //console.log(tpPvalSort[i])
    for(j=0;j<fpPval.length;j++){
      if(parseFloat(fpPval[j])<=parseFloat(tpPvalSort[i])){ //don't forget the parseFloat
        //console.log(fpPval[j])
        fpc++;
      }
    }
    //console.log(fpc);
    //console.log(fpc);
    fdr.push((fpc/possible_fp));
    tdr.push((i+1)/possible_tp); //Since these are sorted, index 0 in the tp (the first one) represents 1/20
  }
 var rates = makeSteps(fdr,tdr);
 return rates;
 //console.log(rates)
 //return [fdr,tdr];
}


function roc(data){
  var fpExp = data[8],
    tpExp = data[9],
    fpPval = data[6],
    tpPval = data[7],//.sort(function(a, b){return a-b}),
    fp05Pval=[],
    fp02Pval=[],
    fp01Pval=[],
    fp005Pval=[],
    fp002Pval=[],
    tp05Pval=[],
    tp02Pval=[],
    tp01Pval=[],
    tp005Pval=[],
    tp002Pval=[];
//console.log(fpExp)
// Get a variable going for each expected frequency there must be a better way to do this but alas
  for(i=0;i<fpExp.length;i++){
    if(parseFloat(fpExp[i])==0.05){
      fp05Pval.push(fpPval[i]);
    }else if(parseFloat(fpExp[i])==0.02){
      fp02Pval.push(fpPval[i]);
    }else if(parseFloat(fpExp[i])==0.01){
      fp01Pval.push(fpPval[i]);
    }else if(parseFloat(fpExp[i])==0.005){
      fp005Pval.push(fpPval[i]);
    }else if(parseFloat(fpExp[i])==0.002){
      fp002Pval.push(fpPval[i]);
    }
  }
  for(i=0;i<tpExp.length;i++){
    if(parseFloat(tpExp[i])==0.05){
      //console.log('found it')
      //console.log(tpPval[i]);
      tp05Pval.push(tpPval[i]);
    }else if(parseFloat(tpExp[i])==0.02){
      tp02Pval.push(tpPval[i]);
    }else if(parseFloat(tpExp[i])==0.01){
      tp01Pval.push(tpPval[i]);
    }else if(parseFloat(tpExp[i])==0.005){
      tp005Pval.push(tpPval[i]);
    }else if(parseFloat(tpExp[i])==0.002){
      tp002Pval.push(tpPval[i]);
    }
  }
//console.log(tp05Pval)
  var roc05 = calcROC(tp05Pval,fp05Pval),
      roc02 = calcROC(tp02Pval,fp02Pval),
      roc01 = calcROC(tp01Pval,fp01Pval),
      roc005 = calcROC(tp005Pval,fp005Pval),
      roc002 = calcROC(tp002Pval,fp002Pval);

return [[roc05[0],roc02[0],roc01[0],roc005[0],roc002[0]],[roc05[1],roc02[1],roc01[1],roc005[1],roc002[1]]]
}




// Functions to make plots

function makeQualPlot(){
var trace1 = {
    x: [30,40],
    y: [10,30],
    mode: 'markers',
    name : "False Positives"

};
var trace2 ={
    x: [36,42],
    y: [36,38],
    mode: 'markers',
    name : "True Positives"
};
var trace3 = {
x: [30,30],
    y: [0,42],
    mode: 'lines',
    name : "Phred cutoff",
    marker : {color : 'black'},
    showlegend: false,
    line: {
      dash: 'dashdot',
      width: 2
    }
};
var trace4 = {
x: [30,42],
    y: [0,0],
    mode: 'lines',
    name : "MapQ cutoff",
    showlegend: false,
    marker : {color : 'black'},
    line: {
      dash: 'dashdot',
      width: 2
    }
};
var data = [trace1,trace2,trace3,trace4];
var layout = {
    title:'Quality Distributions',
    height: 400,
    //width : 350,
    xaxis : {
      title : "Phred",
      ticks: 'outside'
    },
    yaxis : {
      title :"MapQ",
      ticks: 'outside'
    }
  //  width: 350
};
Plotly.newPlot('qualPlot', data, layout);
}
function makeReadPlot(){
  var trace1 = {
      x: [30,40],
      y: [10,30],
      type: 'bar',
      name : "False Positives"
  };
  var trace2 ={
      x: [36,42],
      y: [36,38],
      type: 'bar',
      name : "True Positives"
  };
  var trace3 = {
  x: [0,0],
      y: [0,215],
      mode: 'lines',
      name : "Minimal Read Position",
      showlegend: false,
      marker : {color : 'black'},
      line: {
        dash: 'dashdot',
        width: 2
      }
  };
  var trace4 = {
  x: [125,125],
      y: [0,215],
      mode: 'lines',
      name : "Maximal Read Position",
      showlegend: false,
      marker : {color : 'black'},
      line: {
        dash: 'dashdot',
        width: 2
      }
  };
  var data = [trace1,trace2,trace3,trace4];
  var layout = {
      title:'Distribution on read',
      height: 400,
    //  width: 400,
      barmode : 'group',
      xaxis : {
        title:"Mean positon on read",
        range : [0,125],
        ticks: 'outside'
      }
  };
  Plotly.newPlot('readPlot', data, layout);
}
function makeRocPlot(){
  var trace1 = {
      x: [30,40],
      y: [10,30],
      mode: 'lines',
      name : "5.0%"
  };
  var trace2 ={
      x: [36,42],
      y: [36,38],
      mode: 'lines',
      name : "2.0%"
  };
  var trace3 = {
  x: [0,0],
      y: [0,215],
      mode: 'lines',
      name : "1.0%"
  };
  var trace4 = {
  x: [125,125],
      y: [0,215],
      mode: 'lines',
      name : "0.05%"

  };
  var trace5 = {
  x: [125,125],
      y: [0,215],
      mode: 'lines',
      name : "0.02%"

  };
  var data = [trace1,trace2,trace3,trace4,trace5];
  var layout = {
      title:'ROC curve',
      height: 400,
      //width: 450,
      xaxis : {
        title : "1-Specificity",
        range : [0,0.005],
        zeroline: false,
        ticks: 'outside'
      },
      yaxis : {
        title : "Sensitivity",
        range : [0,1],
        ticks: 'outside'
      }
};
Plotly.newPlot('rocPlot', data, layout);
}
/// Functions to update plots
function upQualPlot(fpPhred,fpMapQ,tpPhred,tpMapQ,phredCutX,phredCutY,mapqCutX,mapqCutY){

  var update = {x : [fpPhred,tpPhred,mapqCutX,phredCutX], y : [fpMapQ,tpMapQ,mapqCutY,phredCutY]};
  //console.log(data)
  Plotly.restyle('qualPlot', update);

}

function upReadPlot(fpPosX,tpPosX,readCutMinX,readCutMaxX,fpPosY,tpPosY,readCutMinY,readCutMaxY){

  var update = {x : [fpPosX,tpPosX,readCutMinX,readCutMaxX], y: [fpPosY,tpPosY,readCutMinY,readCutMaxY]};
  //console.log(data)
  Plotly.restyle('readPlot', update);

}
function upRocPlot(rocData){

  var update = {x : rocData[0], y: rocData[1]}
  Plotly.restyle('rocPlot', update);

}

//make plots to start
makeQualPlot()
makeReadPlot()
makeRocPlot()
function go() {
    var vals = getVals(),
        data = selectData(vals.Titer,vals.MapQ,vals.Phred,vals.ReadPosMin,vals.ReadPosMax, plotData),
        fpPhred = data[0],
        fpMapQ = data[1],
        tpPhred = data[2],
        tpMapQ = data[3],
        fpPosAll = binReads(data[4],4), //5 is the bin size
        tpPosAll = binReads(data[5],4), //5 is the bin size
        fpPosX = fpPosAll.unique(),
        tpPosX = tpPosAll.unique(),
        fpPosY = count(fpPosAll,fpPosX),
        tpPosY = count(tpPosAll,tpPosX)
        mapqCutX = [30,42],
		    phredCutY = [0,42],
		    mapqCutY = [vals.MapQ,vals.MapQ],
		    phredCutX = [vals.Phred,vals.Phred],
        readCutMinY = [0,215],
        readCutMaxY = [0,215],
        readCutMinX = [vals.ReadPosMin,vals.ReadPosMin],
        readCutMaxX = [vals.ReadPosMax,vals.ReadPosMax];
        //console.log(tpPhred);
        //console.log(tpMapQ);
        var rocData=roc(data);
        //console.log(rocData);
    //console.log(fpPhred)
  /// Make qaul plot
        upQualPlot(fpPhred,fpMapQ,tpPhred,tpMapQ,phredCutX,phredCutY,mapqCutX,mapqCutY);
        upReadPlot(fpPosX,tpPosX,readCutMinX,readCutMaxX,fpPosY,tpPosY,readCutMinY,readCutMaxY)
        upRocPlot(rocData)
}
go()
