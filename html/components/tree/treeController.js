app.controller("treeController", function($scope, $http){
	
	$scope.$watch("field.subdata", function(newValue, oldValue){
		
		console.log("[TREE] WATCH", oldValue, newValue);
		
		if(newValue != oldValue){
			
//			console.log("[TREE]", newValue);
			
			var newick = newValue;
			var tree = d3.layout.phylotree().svg(d3.select("#tree")); // .radial(true);
			var phylotree = tree(d3.layout.newick_parser(newick))
			phylotree.layout();
			
//			console.log("[CLOUD] ELEMENT DATA: ", $scope.subdata);
		}
	});
});