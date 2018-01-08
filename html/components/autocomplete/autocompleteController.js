app.controller("autocompleteController", function($scope, $http, $q, $filter, dataService, $timeout){
    
	var self = this;
	
	self.selectedItemChange = function (item) {
		console.log("AUTOCOMPLETE SELECTED", item);
		
		if( $scope.field.data && $scope.field.data.onChange) {
			console.log("GLOBAL CHANGE 4", item);
			
			if($scope.field.data.onChange.field){
				dataService.global[$scope.field.data.onChange.key] = item == undefined ? item : item[$scope.field.data.onChange.field];
			}
			else {
				console.log("GLOBAL CHANGE 5", item);
				dataService.global[$scope.field.data.onChange.key] = item == undefined ? item : item.label || item;
			}
		}
    };
	
	self.querySearch = function(query){
		
		console.log("AUTOCOMPLETE", query, $scope, $scope.field);
		
		if ($scope.field.data && $scope.field.data.url){
			var deferred = $q.defer();
			
			console.log("AUTOCOMPLETE AJAX [GET]", $scope.field, query);
			
			$timeout.cancel(self.filterTextTimeout);

	        self.filterTextTimeout = $timeout(function() {
	        	
	        	$http.get($scope.field.data.url + query).then(
						function(response) {
							
							console.log("AUTOCOMPLETE AJAX [RESULT]", query, response);
							
							if (response.data.length == 1)
								self.selectedItemChange(response.data[0]);
							
							deferred.resolve(response.data);
						},
						function(response) {
						});
	        }, 1000);
			
			return deferred.promise;
		}
		else {
			var results = [];
			
			if(!query) results = $scope.field.data.subdata;
			else if($scope.field.sending) return $scope.field.data.subdata;
			else {
				results = $scope.field.data.subdata.filter(
						function(item)
						{
							return angular.lowercase(item.label).indexOf(angular.lowercase(query)) >= 0;
						});
			}
			
			return results;
		}
	};
	
	if($scope.field.data != undefined){
		
		var keys = [];
		if($scope.field.data.key) keys.push($scope.field.data.key);
		
		for(index in $scope.field.data.templates){
			
			var template = $scope.field.data.templates[index];
			console.log("TEMPLATE2", $scope.field.type, data, template);
			
			if(template.key && keys.indexOf(template.key) == -1) keys.push(template.key);
		}
		if($scope.field.data.onChange && $scope.field.data.onChange.key != undefined && keys.indexOf($scope.field.data.onChange.key) == -1) keys.push($scope.field.data.onChange.key);
		
		for(index in keys)
		{
			var key = keys[index];
			
			if($scope.field.data && $scope.field.data.variable_value)
				$scope.replaceTemplates();
			
//			console.log("ADDING WATCH", $scope.field.type, key, $scope.field, dataService, dataService.global[key], $scope);
//			
//			$scope.$watch(function(){return dataService.global[key];}, function(newValue, oldValue) {
//				console.log("INSIDE WATCH", $scope.field.type, key, newValue, oldValue, $scope);
//				
//			    if (newValue != oldValue){
//			    	
//			    	console.log($scope.field.type, "Variable " + key + " changed from " + oldValue + " to ", newValue, "effective value=", dataService.global[key], " field=", $scope.field);
//			    	
//			    	// $scope.field.data.value = newValue == undefined ? newValue : (newValue.label || newValue);
//			    	// $scope.field.data.value = newValue == undefined ? newValue : newValue[$scope.field.data.onChange.field];
//			    	
//			    	if($scope.field.data && $scope.field.data.onChange != "nothing")
//			    		$scope.update($scope.get_url());
//			    	
//			    	if($scope.field.data.templates)
//			    		$scope.replaceTemplates();
//			    }
//			});
		}
		
		if($scope.field.value) self.querySearch($scope.field.value);
	}
});