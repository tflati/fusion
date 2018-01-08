app.directive("myCustomElement", function() { 
	return { 
		restrict: "E", 
		scope: {
			myElement: "="
		},
		controller : function MyCustomCellController($scope, $mdDialog, $window)
		{
			$scope.onClick = function(ev, elem){
				
				if(elem.action == "window") {
					console.log("Opening a new window with data", elem.card);
					$scope.inputData = elem.card;
					$window.parentScope = $scope;
					var popup = $window.open("/interface-engine/popup", "_blank", "width=800,height=600,left=50,top=50");
				}
				else {
					console.log(ev, elem);
					
					var card = elem.card;
			        
			        $mdDialog.show({
			        	multiple: true,
			        	locals: {data: card},
			            controller: function DialogController($scope, $mdDialog, data) {
			                $scope.row = [data];
			                $scope.closeDialog = function() {
			                  $mdDialog.hide();
			                };
			                
			                $scope.show_dialog = function(ev, card){
			                    
			                    console.log(ev, card);
			                    
			                    $mdDialog.show({
			                    	multiple: true,
			                    	locals: {data: card},
			                        controller: function DialogController($scope, $mdDialog, data) {
			                            $scope.row = [data];
			                            $scope.closeDialog = function() {
			                              $mdDialog.hide();
			                            };
			                          },
			                        templateUrl: 'templates/dialog.html',
			                        parent: angular.element(document.body),
			                        targetEvent: ev,
			                        clickOutsideToClose:true
			                    });
			                };
			              },
			            templateUrl: 'templates/dialog.html',
			            parent: angular.element(document.body),
			            targetEvent: ev,
			            clickOutsideToClose:true
			        });
				}
		    };
		},
		templateUrl: "components/table/my_custom_element/myCustomElement.html"			
	}; 
});