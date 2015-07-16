function createMSA(div_id, url, width){
  var yourDiv = document.getElementById(div_id);
  yourDiv.innerHTML = "<div id='load_msa'>Please wait while we construct the MSA...</div>";

  /* global yourDiv */
  var msa = require("msa");
  var fasta = require("biojs-io-fasta");
  var gff = require("biojs-io-gff");

  var msaDiv = document.createElement('div');
  yourDiv.appendChild(msaDiv);

  var opts = {
    el: msaDiv
  }
  opts.conf = {
    manualRendering: true
  }
  opts.vis = {
    conserv: false,
    overviewbox: false,
    seqlogo: true,
    markers: false,
    leftHeader: false,
    labelName: true,
    labelId: false,
  };
  opts.zoomer = {
    labelNameLength: 160,
    alignmentWidth: width,
    alignmentHeight: 500,
  };

  var m = msa(opts);

  // init msa
  $.ajax({
    url:url,
    dataType: "json",
    success: function(result) {
      m.seqs.reset(result.seqs);
      var features = gff.parseSeqs(result.features);
      m.seqs.addFeatures(features);
      m.render()
      $("#load_msa").html("");
    },
    error: function(jqXHR, textStatus, errorThrown) {
      console.log(textStatus);
      console.log(errorThrown)
    }
  });

  var colorConservation = {}

  // the init function is only called once
  colorConservation.init = function(){
    // you have here access to the conservation or the sequence object
    this.cons = this.opt.conservation();
    this.cons80 = true;
    this.cons50 = true;
  }

  colorConservation.run = function(letter,opts){
    if(this.cons80 && this.cons[opts.pos] > 0.8){
      return "red";
    }
    else if(this.cons50 && this.cons[opts.pos] > 0.5){
      return "blue";
    }
    else{
      return "#000";
    }
  };

  m.g.colorscheme.addDynScheme("colorConservation", colorConservation);
  m.g.colorscheme.set("scheme", "colorConservation");
}