// File Table dierctive
app.directive("myTable", function() {
	return {
		templateUrl : "components/table/myTable.html",
		restrict : "E",
		scope : {
			formData : "="
		},
		controller : function MyDataTableController($scope, $attrs, $timeout, $http, $cookies, dataService) {

			$scope.formData.submit = {
					url : $scope.formData.url,
				};
			if($scope.formData.source) $scope.formData.submit.source = $scope.formData.source;
			
			$scope.formData.fields = [];

			console.log("MYTABLE", $scope.formData, dataService.global);
			
			var fields = [];
			if($scope.formData.submit.source){
				var form = dataService.global[$scope.formData.submit.source];
				for (var el = 0; el < form.elements.length; el++) {
					var element = form.elements[el];
					console.log("ELEMENT", element);
					for (var i = 0; i < element.length; i++) {
						var row = element[i];
						console.log("ROW", row);
						for (var j = 0; j < row.elements.length; j++) {
							var field = row.elements[j];
							if(field.subtype != "form") continue;
							
							console.log("FIELD", field);
							fields.push(field);
						}
					}
				}
			}
			$scope.formData.fields = fields;
			
			console.log("MYTABLE", $scope.formData);
			
//			$scope.formData.results = {};
			$scope.during_call = false;

			$scope.tableUrl = $scope.formData.submit.url;
			if ($scope.tableUrl.slice(-1) != "/") $scope.tableUrl = $scope.tableUrl + "/";

			$scope.filter_visibility = false;

			$scope.showFilter = function() {
				if ($scope.filter_visibility) $scope.filter_visibility = false;
				else $scope.filter_visibility = true;
			};

			$scope.t_property = {
				// 'table-row-id-key': '',
				'column-keys' : []
			};

			$scope.t_id = {};
			$scope.pageSize = 20;
			$scope.filter_list = [];
			
//			for (var i = 0; i < $scope.formData.fields.length; i++) {
//				var field = $scope.formData.fields[i];
//				if(field.key == "submit" && field.results) {
//					$scope.dataReceived(field.results);
//				}
//			}
			
			/* START OF FUNCTIONS */
			$scope.paginatorCallback = function(page, pageSize, options) {

//				$scope.formData.results.hits = [];

				var args = {};
				for (var i = 0; i < $scope.formData.fields.length; i++) {
					var field = $scope.formData.fields[i];
//					console.log("FIELD", field);
					if (!field.data) continue;
					if(!field.key) continue;

					if (angular.isArray(field.data.value)) {
						subargs = []
						for (var j = 0; j < field.data.value.length; j++) {
							var value = field.data.value[j].id;
							if (value == undefined) value = field.data.value[j];
//							console.log("VALUE", value);

							if (value == undefined || value == "undefined" || value == "") value = "ALL";
							subargs.push(value);
						}

						args[field.key] = subargs;
					} else {
						var value = "ALL";
						if (field.data.value && field.data.value.id) value = field.data.value.id;
						else if (field.data.value) value = field.data.value;

						if (field.type == "checkbox") {
							value = field.data.value ? field.data.value : false;
							if (value.value) value = value.value;
						}

						console.log("VALUE", value);

						// if (value == undefined || value ==
						// "undefined" || value == "") value =
						// "ALL";
						args[field.key] = value;
					}
				}

				if (!page) page = 1;
				var offset = (page - 1) * pageSize;

				var postArguments = {
					'offset' : offset,
					'limit' : pageSize,
					'sort' : {
						'field' : 'id',
						'order' : 'ASC'
					},
					'filter' : $scope.filter_list
				};

				for ( var key in args) {
					var value = args[key];
					if (angular.isString(value))
						value = value.replace(/#/g, "_SHARP_");
					postArguments[key] = value;
				}
				console.log("ARGS SENT VIA POST", postArguments);

				$scope.during_call = true;
				dataService.global["num_results"] = "...";
				return $http
						.post($scope.tableUrl, postArguments)
						.then(
								function(result) {

									console.log("TABLE RESULT", result);

									dataService.login();

									$scope.during_call = false;
									var data = result.data;
									data = $scope.dataReceived(data);

									return {
										results : data.hits,
										totalResultCount : result.data.total
									}
								});
			};
			
			$scope.dataReceived = function(data){
				// Base case: no results
				if (!data.structure) data = {
						hits: [],
						totalResultCount : 0,
						structure: {field_list: []}
				};

				// set row table property
				if ($scope.t_property['column-keys'].length == 0) {

					for (var i = 0; i < data.structure.field_list.length; i++) {
						// set column array
						$scope.t_property['column-keys'].push(data.structure.field_list[i].label);

						// set filter array
						if (data.structure.field_list[i].filters.list.length) {
							$scope.filter_list.push({
										"label" : data.structure.field_list[i].label,
										"title" : data.structure.field_list[i].filters.title,
										"filters" : data.structure.field_list[i].filters.list
									});
						}
					}

					$scope.t_column = data.structure.field_list;
				}

				// save the result of the call
				$scope.table_data = data;
//				$scope.formData.results = data;
				
				if($scope.formData.onReceive){
					for(var i=0; i<$scope.formData.onReceive.length; i++){
						var instruction = $scope.formData.onReceive[i];
						if (instruction.key == "num_results" && instruction.action == "write") {
							dataService.global[instruction.key] = data.total;
						}
						else {
							dataService.global[instruction.key] = instruction.value;
						}
					}
				}
				
				return data;
			};

			$scope.getLoadResultsCallback = function(loadPageCallback) {
				$scope.requestAjax = loadPageCallback;
			};

			$scope.applyFilter = function() {
				$scope.requestAjax();
			};

			$scope.resetFilter = function() {
				for (var i = 0; i < $scope.filter_list.length; i++) {
					for (var j = 0; j < $scope.filter_list[i].filters.length; j++) {
						$scope.filter_list[i].filters[j].chosen_value = "";
					}
				}
			};

			$scope.formatFilterLabel = function(label) {
				return label.charAt(0).toUpperCase()
						+ label.slice(1);
			}

//			$scope.getElementById = function(e_id) {
//
//				var true_id = $scope.t_id[e_id];
//
//				// find the true id in the id list
//				for (var i = 0; i < $scope.table_data.hits.length; i++) {
//					if ($scope.table_data.hits[i].id == e_id) {
//						return $scope.table_data.hits[i];
//					}
//				}
//			};
		}
	};
});
