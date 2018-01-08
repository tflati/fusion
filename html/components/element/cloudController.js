app.controller("cloudController", function($scope, $http){
	
//	$scope.init = function(){
//		console.log("[INIT CLOUD]");
		
//		console.log("[CLOUD] SCOPE", $scope);
//		console.log("[CLOUD] PARENT SCOPE", $scope.$parent);
		$scope.width = $scope.$parent.$parent.width || 100;
		$scope.height = $scope.$parent.$parent.height || 100;
		$scope.steps = $scope.$parent.$parent.steps || 7;
	
		$scope.$watch("subdata", function(newValue, oldValue){
			
//			console.log("[CLOUD] WATCH", oldValue, newValue);
			
			if(newValue != oldValue && newValue.hasOwnProperty('items')){
				var aux = [];
				for(var i=0;i<newValue.items.length; i++)
				{
					var item = newValue.items[i];
					aux.push({"text": item[0], "weight": item[1], "html": {title: "Importance: " + item[1]}});
				}
				
				$scope.subdata = aux;
//				console.log("[CLOUD] ELEMENT DATA: ", $scope.subdata);
			}
		});
//	};
});