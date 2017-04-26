app.controller("elementController", function($scope, $http, dataService){
	
	$scope.subdata = [];
	$scope.sending = false;

	$scope.globaldata = dataService.global;
	$scope.$watch(function(){return dataService.global;}, function(newValue) {
	        $scope.globaldata = newValue;
	}, true);
	
	$scope.init = function(data){
		
		$scope.type = data.type;
		$scope.data_source = data.data;
		$scope.numbered = data.numbered;
		
		console.log("ELEMENT METADATA: ", data);
		
		if($scope.data_source.value && $scope.data_source.key)
			dataService.global[$scope.data_source.key] = $scope.data_source.value;
		
		$scope.update($scope.get_url());
		
		if($scope.data_source.key){
			$scope.$watch(function(){return dataService.global[$scope.data_source.key];}, function(newValue, oldValue) {
				
			    if (newValue != oldValue){
			        console.log("Variable " + $scope.data_source.key+ " changed from " + oldValue + " to " + newValue + " (effective value="+dataService.global[$scope.data_source.key]+")");
			        $scope.update($scope.get_url());
			    }
			});
		}
	};
	
	$scope.get_url = function(){
		var template = $scope.data_source.template;
		
		var value = dataService.global[$scope.data_source.key];		
		var url = $scope.data_source.url;
		if(!template) return url;
		else return url.replace(template, value);
	};
	
	$scope.update = function(url){
		
		if(!url.endsWith("/")) url += "/";
		
		console.log("AJAX TO", url);
		$scope.sending = true;
		$http.get(url).then(
			function(response)
    		{
				$scope.sending = false;
				console.log("SUCCESS IN GETTING DATA FROM " + url, response.data);
				$scope.subdata = response.data.details;
				
				console.log("ELEMENT DATA: ", $scope.subdata);
			},
			function myError(response) {
				$scope.sending = false;
            	console.log("ERROR IN GETTING DATA FROM " + url, response);
			}
		);
	};
	
	$scope.onClick = function(evt){
		
		var value = $scope.subdata.header[evt[0]._index];
		console.log("CLICK EVENT", evt, value);
//	    var chr_for_query = chr.replace("chr", "").trim();
//	    self.clicked_chromosome = chr_for_query;
//	    self.load_subdata(self.info.links.statistics_single_chromosome + chr_for_query + "/");
		
		var listener = $scope.data_source.onClick;
		if(listener.action == "write") {
			
		$scope.$apply(function() {
			console.log("GLOBAL:", dataService);
			dataService.global[listener.key] = value;
			console.log("CHANGING VALUE OF VARIABLE '" + listener.key + "' TO " + value + " (real value: "+dataService.global[listener.key]+")");
		    });
		}
	};
});