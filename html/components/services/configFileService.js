app.service("configFileService", function($http, $q, messageService, dataService, Analytics){
	
	var self = this;
	
	self.info = undefined;
	self.ajax2forms = {};
	
	self.get_config_file = function(){
		
		if(self.info == undefined) {
			console.log("self.info is undefined. Loading...");
			return self.load_config_file();
		}
		
		console.log("The config file has already been downloaded.");
		
		var deferred = $q.defer();
		deferred.resolve(self.info);		
		return deferred.promise;
	};
	
	self.load_config_file = function(){
		
		console.log("MAKING AJAX FOR CONFIG FILE");
		
		return $http.get('config.json').then(function(response) {
			
			console.log("File di configurazione caricato correttamente.");
			// messageService.showMessage('File di configurazione caricato correttamente.');
			console.log("FILE CONFIG OK");
			
	        self.info = response.data;
	        
	        console.log(Analytics);
	        Analytics.configuration.accounts[0].tracker = self.info.analytics_tracking_ID;
	        Analytics.registerScriptTags();
	        Analytics.registerTrackers();
	        
	        self.info.image.percentage_width_original = self.info.image.percentage_width;
	        console.log(self.info);
	        
	        if(self.info.forms)
		        for(var f=0; f<self.info.forms.length; f++)
		        {
		        	var form = self.info.forms[f];
		        	dataService.global[form.pageID] = form;
		        }
	        
	        if(self.info.pages)
		        for(var f=0; f<self.info.pages.length; f++)
		        {
		        	var page = self.info.pages[f];
		        	dataService.global[page.pageID] = page;
		        }
	        
	        return self.info;
	        
		}, function(response){
			console.log(response, 'Impossibile trovare il file di configurazione. (message: "' + response.statusText + '", code: '+response.status+')');
			// messageService.showMessage('Impossibile trovare il file di configurazione. (message: "' + response.statusText + '", code: '+response.status+')', "error", "Error");
			
			return undefined;
		});
	};
});