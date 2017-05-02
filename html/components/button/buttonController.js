app.controller("buttonController", function($scope){
	
	$scope.init = function(data){
            
            console.log("BUTTON CREATION", data);
            
		$scope.type = data.type;
		$scope.label = data.label;
		$scope.icon = data.icon;
		$scope.size = data.size;
		$scope.url = data.url;
	};
});