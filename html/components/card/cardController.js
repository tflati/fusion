app.controller("cardController", function($scope, $window, $rootScope, dataService){
	
	var self = this;
	
	$scope.sending = false;
	
	$scope.init = function(data){
		
		console.log("CARD CONTROLLER", data.title, data);
		
//		console.log("CARD SCOPE", data, $scope);
//		console.log("CARD DATA", data);
		$scope.card = data;
		
//		$scope.title = data.title;
//		$scope.elements = data.elements;
//		$scope.subtitle = data.subtitle;
//		$scope.type = data.type;
//		$scope.alignment = data.alignment;
//		$scope.width = data.width;
//		$scope.compact = data.compact,
//		$scope.footer = data.footer;
//		$scope.elements = data.elements;
//		$scope.layout = data.layout;
		
		if(!$scope.card.show) $scope.card.show = true;
		
//		$scope.show = true;
//		if (data.show == false) $scope.show = false;
		
//		if (data.show == "after_send"){
//			$scope.card.show = $rootScope.search_started;
//
//			$scope.$watch(function(){return $rootScope.search_started;}, function(newValue, oldValue) {
//				console.log("SEARCH STARTED WATCH", $rootScope.search_started);
//				$scope.card.show = newValue;
//			}, true);
//		}
//		
//		if (data.show == "before_send"){
//			$scope.card.show = !$rootScope.search_started;
//
//			$scope.$watch(function(){return $rootScope.search_started;}, function(newValue, oldValue) {
//				console.log("SEARCH STARTED WATCH", $rootScope.search_started);
//				$scope.card.show = !newValue;
//			}, true);
//		}
		
		for(var i=0; i < $scope.card.show.length; i++)
		{
			var condition = $scope.card.show[i];
			var arg1 = condition.arg1;
			if(arg1.key){
				console.log("CARD SHOW ADDING WATCH TO ", arg1);
				$scope.$watch(function(){return dataService.global[arg1.key];}, function(newValue, oldValue) {
					console.log("WATCH CARD SHOW", oldValue, newValue);
					condition.arg1.value = newValue;
//					if(newValue != oldValue)
						$scope.updateShow(condition);
				});
			}
			
			var op = condition.op;
			var arg2 = condition.arg2;
			if(arg2.key){
				console.log("CARD SHOW ADDING WATCH TO ", arg2);
				$scope.$watch(function(){return dataService.global[arg2.key];}, function(newValue, oldValue) {
					console.log("WATCH CARD SHOW", oldValue, newValue);
					condition.arg2.value = newValue;
//					if(newValue != oldValue)
						$scope.updateShow(condition);
				});
			}
			
//			$scope.updateShow(condition);
		}
		
		console.log("CARD DATA: ", data);
	};
	
	$scope.updateShow = function(condition){
		console.log("WATCH CARD SHOW UPDATE", condition);
		
		var arg1 = condition.arg1;
		if(arg1 == "undefined") arg1 = undefined;
		if(condition.arg1.fx == "value") arg1 = arg1.value;
		else if(condition.arg1.fx == "length") arg1 = length(arg1.value);
		
		var arg2 = condition.arg2;
		if(arg2 == "undefined") arg2 = undefined;
		if(condition.arg2.fx == "value") arg2 = arg2.value;
		else if(condition.arg2.fx == "length") arg2 = length(arg2.value);
		
		var show = false;
		if(condition.op == ">=") show = arg1 >= arg2;
		else if(condition.op == ">") show = arg1 > arg2;
		else if(condition.op == "<") show = arg1 < arg2;
		else if(condition.op == "<=") show = arg1 <= arg2;
		else if(condition.op == "==") show = arg1 == arg2;
		
		$scope.card.show = show;
	};
	
	$scope.footerFunction = function(fx){
		if(fx == 'total') return ($scope.elements ? $scope.elements[0].data.length : 0);
		else return "UNKNOWN";
	};
	
	if($window.opener && $window.opener.parentScope && $window.opener.parentScope.inputData)
		$scope.init($window.opener.parentScope.inputData);
});