app.controller("buttonController", function($scope){
	
	$scope.init = function(data){
		$scope.type = data.type;
		$scope.label = data.label;
		$scope.icon = data.icon;
		$scope.size = data.size;
		$scope.url = data.url;
		
		console.log("BUTTON CREATION", data);
	};
});