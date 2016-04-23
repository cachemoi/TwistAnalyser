/**
 * Created by Maximilien Rothier Bautzer
 */
function Unit (queryName, specie, strand1, strand2, strand3) {

    this.qname = queryName;
    this.specie = specie;
    this.strand1 = strand1;
    this.strand2 = strand2;
    this.strand3 = strand3
}

function Options(UP, PDB, UI, DOM, GXY_num, inc, SP, OV, ND){
    this.UP = UP;
    this.PDB =PDB;
    this.UI = UI;
    this.ND = ND;
    this.DOM = DOM;
    this.GXY_num = GXY_num;
    this.inc = inc;
    this.SP = SP;
    this.OV = OV
}

function FileInfo(){
    this.filePath = "http://localhost/TwistAnalyser/";
    this.filename = "";

    this.setFile = function (filename, file){
        this.filePath = "http://localhost/TwistAnalyser/" + file +"/" +filename;
        this.filename = filename;
    };

    this.reset = function (){
        this.filePath = "http://localhost/TwistAnalyser/";
        this.filename = "";
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

function checkParameters() {

    options = new Options();

    if ($('#inc').is(":checked")){
        options.inc = true
    }
    if ($('#SP').is(":checked")){
        options.SP = true
    }
    if ($('#OV').is(":checked")){
        options.OV = true
    }

    if ($('#GXY_num_check').is(":checked")){
        GXY_num = $('#GXY_num').val();
        options.GXY_num = GXY_num;
    }

    alert("Data sent")
}

function userInputFile() {

    var isUP = $('#UPfileType').is(":checked");
    var isPDB = $('#PDBfileType').is(":checked");

    unit = new Unit();

    unit.qname = $('#qnamef').val();

    unit.specie = $('#FileSpecie').val();

    if (isUP){

        var i = fileInputResults.indexOf('\n');

        asArray = [fileInputResults.slice(0,i), fileInputResults.slice(i+1)];

        deformat = asArray[1].replace(/(\r\n|\n|\r)/gm,"");

        unit.strand1 = deformat;
        unit.strand2 = deformat;
        unit.strand3 = deformat;

    }

    if (isPDB){

        fileInputResults.replace(/(\r\n|\n|\r)/gm,"");

        var count = (fileInputResults.match(/SEQUENCE/g) || []).length;

        asArray = fileInputResults.split(/\n/);

        switch(count){
            case 1:
                unit.strand1 = asArray[1];
                unit.strand2 = asArray[1];
                unit.strand3 = asArray[1];
                break;
            case 2:
                unit.strand1 = asArray[1];
                unit.strand2 = asArray[3];
                unit.strand3 = asArray[3];
                break;
            case 3:
                unit.strand1 = asArray[1];
                unit.strand2 = asArray[3];
                unit.strand3 = asArray[5];
                break;
        }
    }

    if (!(isPDB || isUP)) {
        alert("Please choose your file type");
    }
    if(fileInputResults==null){
        alert("Please choose a file");
    }

    options.UI = unit.qname + "?" + unit.specie + "?" + unit.strand1 + "|" + unit.strand2 + "|" + unit.strand3;

    console.log(options.UI)
}

function userInputText() {

    // the error trapping here needs way more work


    var specie = $("#Specie").val();
    specie.replace(/" "/g,'').toUpperCase().replace(/(\r\n|\n|\r)/gm,"");
    var strand1 = $("#s1").val();
    strand1.replace(/" "/g,'').toUpperCase().replace(/(\r\n|\n|\r)/gm,"");
    var strand2 = $("#s2").val();
    strand2.replace(/" "/g,'').toUpperCase().replace(/(\r\n|\n|\r)/gm,"");
    var strand3 = $("#s3").val();
    strand3.replace(/" "/g,'').toUpperCase().replace(/(\r\n|\n|\r)/gm,"");

    if (specie == null || strand1== null || strand2== null || strand3== null){
        console.log("loop running");
        alert("please fill all input boxes")
    }

    unit = new Unit(specie, strand1, strand2, strand3);

    if (specie == null || strand1== null || strand2== null || strand3== null){
        alert("please fill all input boxes")
    }

    options.qname = $("#qname").val();

    options.UI = (options.qname+ "?" + specie + "?" +  strand1 + "|" +  strand2 + "|" +  strand3);

}

function processPDBID(){

    options.PDB = true;

    if ($('#DOMPDB').is(":checked")){
        options.DOM = true;
    }

    ID=$("#PDB").val();

    options.ND = ID.split(" ");

}

function processUPID(){

    options.UP = true;

    if ($('#DOMUP').is(":checked")){
        options.DOM = true;
    }

    ID=$("#UP").val();

    options.ND = ID.split(" ");

}

function serverComs(){

    //the codes need to go in ND, but it only works with uniprot for now
    console.log(JSON.stringify(options));
    return $.ajax({
        url: "http://localhost/cgi-bin/python.py",
        type: "POST",
        data: JSON.stringify(options),
        dataType: "text",
        success: function(response){
            console.log("it's alive");
            handleResponse (response)
        },
        error: function(){
            alert('request failed');
        }
    })
}

function handleResponse (response) {
    console.log("response recieved");
    console.log("the response printed");
    console.log(response);
    console.log("the response as json");
    console.log(typeof response);
    commandOUT = JSON.parse(response);
    console.log(commandOUT);
    console.log("the codes are");
    console.log(commandOUT.P);

    var fileinfo = new FileInfo();

    resultDisplay = $("#results");

    var PLength = commandOUT.P.length;
    for (var i = 0; i < PLength; i++) {

        fileinfo.setFile(commandOUT.P[i] + ".txt", "Data");
        resultDisplay.prepend("<br><a href=\"" + fileinfo.filePath + "\" download=\""+ commandOUT.P[i] + "\">Download results</a>");
        fileinfo.reset();

        fileinfo.setFile(commandOUT.P[i] + ".png", "Graphs");
        filepath = "<img src=\"" + fileinfo.filePath + " \" alt=\""+ commandOUT.P[i] + "\">";
        resultDisplay.prepend("<img src=\"" + fileinfo.filePath + " \" alt=\""+ commandOUT.P[i] + "\">");
        resultDisplay.prepend("<h2>" + commandOUT.P[i] + "</h2>");
        fileinfo.reset();
    }
    resultDisplay.prepend("<h1> Individual proteins results</h1>");

    if (commandOUT.SPMD){
        var SPLength = commandOUT.SP.length;
        for (var y = 0; y < SPLength; y++) {
            console.log(commandOUT.SP);
            fileinfo.setFile(commandOUT.SP[y] + ".txt", "Data");
            resultDisplay.prepend("<br><a href=\"" + fileinfo.filePath + "\" download=\""+ commandOUT.SP[y] + "\">Download results</a>");
            fileinfo.reset();

            fileinfo.setFile(commandOUT.SP[y] + " Metadata.png", "Graphs");
            resultDisplay.prepend("<img src=\"" + fileinfo.filePath + "\" alt=\""+ commandOUT.SP[y] + "\">");
            resultDisplay.prepend("<h2>" +  commandOUT.SP[y] + "</h2>");
            fileinfo.reset();
        }

        resultDisplay.prepend("<h1> Individual species results</h1>");

        fileinfo.setFile("Species Metadata.png", "Graphs");
        console.log(fileinfo.filePath);
        resultDisplay.prepend("<img src=\"" + fileinfo.filePath + "\" alt=\"SPM\">");
        resultDisplay.prepend("<h1> Species Metadata</h1>");
        fileinfo.reset();
    }

    if (commandOUT.MD){

        fileinfo.setFile("Overall.txt", "Data");
        console.log(fileinfo.filePath);
        resultDisplay.prepend("<br><a href=\"" + fileinfo.filePath + "\" download=\"Overall metadata\">Download results</a>");
        fileinfo.reset();

        fileinfo.setFile("Metadata.png", "Graphs");
        resultDisplay.prepend("<img src=\"" + fileinfo.filePath + " \" alt=\"Meta\">");
        resultDisplay.prepend("<h1>Metadata</h1>");
        fileinfo.reset();
    }

    resultDisplay.prepend("<h3>Overall average angle: " + commandOUT.OAA.toFixed(2) + "</h3>");

    resultDisplay.prepend("<h3>Overall average translation: " + commandOUT.OAT.toFixed(2) + "</h3>");

    showResults();
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
    serverComs().done(handleResponse(response));
}