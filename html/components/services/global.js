app.filter('pretty', function() {

  return function(input) {
	  return input.replace(/\n/g, "<br/>");
  }
});

app.service("dataService", function($cookies) {

    this.global = {};
    
    this.global["username"] = $cookies.get('username');
    this.username = $cookies.get('username');
    this.loggedIn = $cookies.get('logged_in');
    this.loginToken = $cookies.get('login_token');
    
    this.login = function(){
    	this.global["username"] = $cookies.get('username');
        this.username = $cookies.get('username');
        this.loggedIn = $cookies.get('logged_in');
        this.loginToken = $cookies.get('login_token');
        
//    	dataService.username = $cookies.get('username');
//		dataService.global["username"] = $cookies.get('username');
//		dataService.loggedIn = $cookies.get('logged_in');
//		dataService.loggedToken = $cookies.get('login_token');
		console.log("COOKIE UPDATE", this.loggedIn, this.loggedToken);    	
    }
});