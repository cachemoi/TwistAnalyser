/**
 * Created by Maximilien Rothier Bautzer
 */

function Strand(input) {

    this.aminoSequence= input;
    this.strandPLocation = [];

    this.splitStrand1 = function () {
        this.aminoSequence = this.aminoSequence.slice(2);
    };

    this.splitStrand2 = function () {
        this.aminoSequence = this.aminoSequence.slice(1);
        this.aminoSequence = this.aminoSequence.slice(0, -1)
    };

    this.splitStrand3 = function(){
        this.aminoSequence = this.aminoSequence.slice(0, -2);
    };

    this.lookForP = function (){
        for ( var i = 0; i < this.aminoSequence.length; i++ ) {
            if (this.aminoSequence.charAt(i)== "P"){
                this.strandPLocation.push(1);
            }else{
                this.strandPLocation.push(0);
            }
        }
    }
}

// merging all the processed strands

function allThreeStrands(strand1, strand2, strand3){

    this.totalPValues =[];

    this.strand1 = strand1;
    this.strand2 = strand2;
    this.strand3 = strand3;

    for(var i = 0; i < this.strand1.length; i++){
            this.totalPValues.push(this.strand1[i] + this.strand2[i] + this.strand3[i]);
    }
}

window.onload = function() {
    fileInput = document.getElementById('inputFile');
    fileDisplayArea = document.getElementById('fileDisplayArea');

    fileInput.addEventListener('change', function(e) {
        var file = fileInput.files[0];
        var textType = /text.*/;

        if (file.type.match(textType)) {
            var reader = new FileReader();

            reader.onload = function(e) {
                fileDisplayArea.innerText = reader.result;
                fileInputResults = reader.result;
            };
            reader.readAsText(file);
        } else {
            fileDisplayArea.innerText = "File not supported!"
        }
    });
};

function userInputFile() {

    strandUser1 = new Strand ("");
    strandUser2 = new Strand ("");
    strandUser3 = new Strand ("");

    lines = fileInputResults.split('\n');
    var count = 1;
    for (var i = 0; i<lines.length; i++){
        switch (count){
            case 1:
                strandUser1.aminoSequence += (lines[i]).trim();
                count++;
                console.log("case1");
                break;
            case 2:
                strandUser2.aminoSequence += (lines[i]).trim();
                count++;
                console.log("case2");
                break;
            case 3:
                strandUser3.aminoSequence += (lines[i]).trim();
                i++;
                console.log("case3");
                count =1;
                break;
            default:
                console.log("case4");
                break;
        }
    }

    strandUser1.splitStrand1();
    strandUser1.lookForP();

    strandUser2.splitStrand2();
    strandUser2.lookForP();


    strandUser3.lookForP();
    strandUser3.splitStrand3();

    all3UserStrands = new allThreeStrands(strandUser1.strandPLocation,strandUser2.strandPLocation,strandUser3.strandPLocation);

}

function userInput() {
    strandUser1 = new Strand (document.getElementById("s1").value);
    strandUser1.splitStrand1();
    strandUser1.lookForP();

    strandUser2 = new Strand (document.getElementById("s2").value);
    strandUser2.splitStrand2();
    strandUser2.lookForP();

    strandUser3 = new Strand (document.getElementById("s3").value);
    strandUser3.lookForP();
    strandUser3.splitStrand3();

    all3UserStrands = new allThreeStrands(strandUser1.strandPLocation,strandUser2.strandPLocation,strandUser3.strandPLocation);

}

function stepDistribCalc() {
    angles =[];
    translation = [];
    userPValues = all3UserStrands.totalPValues;

    freqTwoTwo = 0;
    freqTwoOne = 0;
    freqTwoNull = 0;
    freqOneOne = 0;
    freqOneNull = 0;
    freqNullNull = 0;

    for(var i = 0; i < (userPValues.length)-1; i++){

        //TYPE 1 collagen

        switch ("" + userPValues[i] + userPValues[i+1]){

            case "22":
                angles.push(-102.6);
                translation.push(2.84);
                freqTwoTwo++;
                break;
            case "21":
                angles.push(-104.5);
                translation.push(2.86);
                freqTwoOne++;
                break;
            case "12":
                angles.push(-104.5);
                translation.push(2.86);
                freqTwoOne++;
                break;
            case "20":
                angles.push(-105.0);
                translation.push(2.86);
                freqTwoNull++;
                break;
            case "02":
                angles.push(-105.0);
                translation.push(2.86);
                freqTwoNull++;
                break;
            case "11":
                angles.push(-105.0);
                translation.push(2.86);
                freqOneOne++;
                break;
            case "10":
                angles.push(-105.9);
                translation.push(2.89);
                freqOneNull++;
                break;
            case "01":
                angles.push(-105.9);
                translation.push(2.89);
                freqOneNull++;
                break;
            case "00":
                angles.push(-108.1);
                translation.push(2.90);
                freqNullNull++;
                break;
            default:
                angles.push("n");
                translation.push("n");
                break;
            }
    }
    freqTotal = freqTwoTwo + freqTwoOne + freqTwoNull + freqOneOne + freqOneNull + freqNullNull;
}

function calculateAverages(){
    averageAngle = 0;
    averageTranslation = 0;
    totalAngles = 0;
    totalTranslation = 0;
    for(var i = 0; i < angles.length; i++){
        totalAngles+=angles[i]
    }

    averageAngle =totalAngles/angles.length;

    for(var x = 0; x < translation.length; x++){
        totalTranslation+=translation[x];
    }

    averageTranslation = totalTranslation/translation.length;

}

function publishTotalPValues(){
    document.getElementById("PDistrib").innerHTML= all3UserStrands.totalPValues.toString();
}
function publishDegreeDistrib(){

    document.getElementById("results").style.display="block";
    document.getElementById("AA").innerHTML = averageAngle.toFixed(2);
    document.getElementById("AT").innerHTML = averageTranslation.toFixed(2);

    console.log(strandUser1.aminoSequence);
    console.log(strandUser2.aminoSequence);
    console.log(strandUser3.aminoSequence);

    console.log(strandUser1.aminoSequence.length);
    console.log(strandUser2.aminoSequence.length);
    console.log(strandUser3.aminoSequence.length);

    document.getElementById("freqTwoTwo").innerHTML = freqTwoTwo;
    document.getElementById("freqTwoOne").innerHTML = freqTwoOne;
    document.getElementById("freqTwoNull").innerHTML = freqTwoNull;
    document.getElementById("freqOneOne").innerHTML = freqOneOne;
    document.getElementById("freqOneNull").innerHTML = freqOneNull;
    document.getElementById("freqNullNull").innerHTML = freqNullNull;
    document.getElementById("freqTotal").innerHTML = freqTotal;

}
function drawAngleChart() {

    // Create the data table.
    var data = new google.visualization.DataTable();
    data.addColumn('number', 'Position');
    data.addColumn('number', 'angle');
    for(var i = 0; i < angles.length; i++){
        data.addRows([
            [i+1, parseFloat(angles[i].toFixed(2))],
        ]);
    }

    var options = {
        titlePosition: 'none',
        legend:'none',
        lineWidth: 1,
        explorer: true,
        hAxis: {
            titleTextStyle: {
                fontSize: 20,
                bold: true,
                italic: false
            },
            title: 'Position',
            fontSize: 15,
            textStyle: {
                italic: false
            }
        },
        vAxis: {
            title: 'K(°)',
            titleTextStyle: {
                fontSize: 20,
                bold: true,
                italic: false
            },
            textStyle: {
                fontSize: 10,
                italic: false
            }
        }
    };

    var chart = new google.visualization.LineChart(document.getElementById('angleChart'));
    chart.draw(data, options);

}
function drawTranslationChart() {

    // Create the data table.
    var data = new google.visualization.DataTable();
    data.addColumn('number', 'Position');
    data.addColumn('number', 'translation');
    for(var i = 0; i < translation.length; i++){
        data.addRows([
            [i+1, parseFloat(translation[i].toFixed(2))],
        ]);
    }

    var options = {
        titlePosition: 'none',
        legend:'none',
        lineWidth: 1,
        explorer: true,
        hAxis: {
            title: 'Position',
            titleTextStyle: {
                fontSize: 15,
                fontName: 'Arial',
                bold: true,
                italic: false
            },
            textStyle: {
                fontSize: 10,
                italic: false
            }
        },
        vAxis: {
            title: 't (‎Å)',
            titleTextStyle: {
                fontSize: 15,
                fontName: 'Arial',
                bold: true,
                italic: false
            },
            textStyle: {
                fontSize: 10,
                italic: false
            }
        }
    };

    var chart = new google.visualization.LineChart(document.getElementById('translationChart'));
    chart.draw(data, options);

}
function serverComs(){
    return $.ajax({
        url: "http://localhost/cgi-bin/python.cgi",
        type: "POST",
        data: JSON.stringify({"UI": false,
            "PDB": [],
            "UP": [],
            "Fibril domains": false,
            "strand2": "GXYGXYGXYGXYGXYGXYGXYGXY",
            "strand3": "GXYGXYGXYGXYGXYGXYGXYGXY",
            "specie" : "Max Bautzer"
        }),
        dataType: "json"
    })
}
        
function handleResponse (response) {
    alert(response);
    var answer = JSON.parse(response);
    fileDisplayArea.innerText = answer.specie;
    console.log("test");
    console.log(response);

    showResults()
}

function showResults(){
    $("#results").show();
    $("#input").hide()
}
function showInput(){
    $("#results").hide();
    $("#input").show();
}

function sendInput(){
    serverComs().done(handleResponse);
}