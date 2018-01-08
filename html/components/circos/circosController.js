app.controller("circosController", function($scope, $http){
	
	var self = this;
	
	console.log("CIRCOS", "init", $scope);
	
	self.circos = { value: "CCLE_001", values: [], events: [ "LINK01", {
		  LinkRadius: 60,
		  LinkFillColor: "#F26223",
		  LinkWidth: 3,
		  displayLinkAxis: true,
		  LinkAxisColor: "#B8B8B8",
		  LinkAxisWidth: 0.5,
		  LinkAxisPad: 3,
		  displayLinkLabel: true,
		  LinkLabelColor: "red",
		  LinkLabelSize: 13,
		  LinkLabelPad: 8,
		}],
		genome: [
		    ["1" , 249250621],
		    ["2" , 243199373],
		    ["3" , 198022430],
		    ["4" , 191154276],
		    ["5" , 180915260],
		    ["6" , 171115067],
		    ["7" , 159138663],
		    ["8" , 146364022],
		    ["9" , 141213431],
		    ["10" , 135534747],
		    ["11" , 135006516],
		    ["12" , 133851895],
		    ["13" , 115169878],
		    ["14" , 107349540],
		    ["15" , 102531392],
		    ["16" , 90354753],
		    ["17" , 81195210],
		    ["18" , 78077248],
		    ["19" , 59128983],
		    ["20" , 63025520],
		    ["21" , 48129895],
		    ["22" , 51304566],
		    ["X" , 155270560],
		    ["Y" , 59373566]
	 ]};
	
	$scope.$watch("field.subdata", function(newValue, oldValue){
		if(newValue != oldValue && newValue.hasOwnProperty('items')){
			
			$("#biocircos").empty();
			self.circos.events[2] = [];
			
			var circos_events = [];
	        
	        for(var i=0; i<newValue.items.length; i++)
	        {
	                var item = newValue.items[i];
	                
	                // Prototype: {fusion: "FGFR3--TACC3", g1chr: 4, g1start: 1795662, g1end: 1808986, g1name: "FGFR3", g2chr: 4, g2start: 1723217, g2end: 1746905, g2name: "TACC3"},
	                
	                if(i<10) console.log(item);
	                
	                if(!$scope.field.data.hide_names){
		                item.fusion = item.g1 + "--" + item.g2;
		                item.g1name = item.g1;
		                item.g2name = item.g2;
	                }
	                delete item.g1
	                delete item.g2
	                
	                item.g1start = parseInt(item.g1start);
	                item.g1end = item.g1start
	                item.g2start = parseInt(item.g2start);
	                item.g2end = item.g2start
	                item.g1chr = item.g1chr.replace("chr", "");
	                item.g2chr = item.g2chr.replace("chr", "");
	                
	                circos_events.push(item);
	        }                              
	        
	        console.log("CIRCOS EVENTS", circos_events);
	        circos_events = circos_events.slice(0, 50);
	        
	        self.circos.events[2] = circos_events;
	        self.circos.events[1].LinkRadius = parseInt($scope.field.data.radius) || 60;
	        
	        BioCircos01 = new BioCircos(self.circos.events, self.circos.genome, {
	            target : "biocircos",
	            svgWidth : parseInt($scope.field.data.width) || 600,
	            svgHeight : parseInt($scope.field.data.height) || 400,
	            innerRadius: parseInt($scope.field.data.inner_radius) || 160,
	            outerRadius: parseInt($scope.field.data.outer_radius) ||180,
	            genomeFillColor: ["#FFFFCC", "#CCFFFF", "#FFCCCC", "#CCCC99","#0099CC", "#996699", "#336699", "#FFCC33","#66CC00"],
	            LINKMouseEvent : true,
	            LINKMouseClickDisplay : true,
	            LINKMouseClickOpacity : 1.0,
	            LINKMouseClickStrokeColor : "red",
	            LINKMouseClickStrokeWidth : 6,
	            LINKLabelDragEvent : false,
	            });
	            BioCircos01.draw_genome(BioCircos01.genomeLength);
			
			self.circos.values = newValue;
		}
	});
});