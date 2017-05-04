app.controller("controller", function($http, $scope, $mdDialog, $timeout, $mdSidenav){

                self = this;
                
                self.page = 'templates/main.html';
                self.info = {};
                self.header = {show_logos: true};
                
                self.sending = true;
                
                self.form;
                self.form_results = [];
                self.sending = false;

                self.circos = { value: "CCLE_001", values: [], events: [ "LINK01", {
                            LinkRadius: 60,
                            LinkFillColor: "#F26223",
                            LinkWidth: 3,
                            displayLinkAxis: true,
                            LinkAxisColor: "#B8B8B8",
                            LinkAxisWidth: 0.5,
                            LinkAxisPad: 3,
                            displayLinkLabel: true,
                            LinkLabelColor: "red",
                            LinkLabelSize: 13,
                            LinkLabelPad: 8,
                        }],
                        genome: [
                            ["1" , 249250621],
                            ["2" , 243199373],
                            ["3" , 198022430],
                            ["4" , 191154276],
                            ["5" , 180915260],
                            ["6" , 171115067],
                            ["7" , 159138663],
                            ["8" , 146364022],
                            ["9" , 141213431],
                            ["10" , 135534747],
                            ["11" , 135006516],
                            ["12" , 133851895],
                            ["13" , 115169878],
                            ["14" , 107349540],
                            ["15" , 102531392],
                            ["16" , 90354753],
                            ["17" , 81195210],
                            ["18" , 78077248],
                            ["19" , 59128983],
                            ["20" , 63025520],
                            ["21" , 48129895],
                            ["22" , 51304566],
                            ["X" , 155270560],
                            ["Y" , 59373566]
                    ]};
                
                self.print = function(){
                        console.log(self.info.forms);
                };
                
                self.load_form = function(group, option){
                        
                        // Remove previous form inputs...
                        for(var f=0; f<self.form.fields.length;)
                        {
                                //console.log("FORM FIELD: ", f, self.form.fields[f]);
                                
                                if(self.form.fields[f].parent_group_id == group.group_id)
                                {
                                        //console.log("Removing form field " + f, self.form.fields[f]);
                                        // toRemove.push(f);
                                        self.form.fields.splice(f, 1);
                                }
                                else f++;
                        }
                        
                        // ... and add newly selected form inputs
                        for(var f=0; f<option.form.fields.length; f++)
                        {
                                // console.log("Adding form field " + f, option.form.fields[f]);
                                option.form.fields[f].parent_group_id = group.group_id;
                                var newForm = option.form.fields[f];
                                
                                self.form.fields.push(newForm);
                                
                                /*
                                if(newForm.type === "checkbox" && newForm.value === "true")
        {
                                        newForm.checked = true;
        }
        else if(newForm.type === "select")
        {
                selectIndex = f;
                
                console.log("GETTING VALUES FOR SELECT: " + newForm.values);
                
                $http.get(newForm.values).then(function(response) {
                        console.log("SELECT VALUES:",response);
                        console.log("SELECT INDEX:", selectIndex);
                        
                        var selectField = option.form.fields[selectIndex];
                        selectField.values = response.data;
                        console.log("NEW VALUES A: ", selectField);
                        console.log("NEW VALUES B: ", option.form.fields[selectIndex]);
                }, function(response){console.log("COULD NOT GET VALUES FOR SELECT FIELD", response);});
        }
                                */
                        }
                };
                
                self.goTo = function(item){
                        url = item.url
                        console.log("Want to go to " + url);
                        
                        // Check if the url is a form name
                        var isForm = false;
                        for(var f=0; f<self.info.forms.length; f++)
                        {
                                var form = self.info.forms[f];
                                if(form.name == url)
                                {
                                        isForm = true;
                                        self.form = form;
                                        
                                        console.log(self.form);
                                        
                                        break;
                                }
                        }
                        
                        if(isForm) {
                                self.form_results = [];
                                self.page = 'templates/form.html';
                                self.info.image.percentage_width = self.info.image.percentage_width_original / 2;
                                self.header.show_logos = false;
                        }
                        else {
                                if(url == "home") url = "main";
                                
                                self.page = 'templates/'+url+'.html';
                                
                                if(item.action) {
                                        self.info.image.percentage_width = self.info.image.percentage_width_original / 2;
                                        self.header.show_logos = false;
                                        self.load_data(item.action);

                                        console.log("CONTROLLER", self.info);
                                }
                                
                                if(url == "main") {
                                        self.page = 'templates/main.html'; // window.location = url;
                                        self.info.image.percentage_width = self.info.image.percentage_width_original;
                                        self.header.show_logos = true;
                                }
                        }
                };
                
                self.show_circos = function(url){
                        $("#biocircos").empty();
                        self.sending = true;
                        self.circos.events[2] = [];
                        
                        $http.get(url).then(function(response) {
                                self.sending = false;
                                
                                console.log("SHOWING CIRCOS DATA: SUCCESS from url", url, response.data.rows);
                                
                                var event = {fusion: "ABL1--BCR", g1chr: 9, g1start: 133589707, g1end: 133589707, g1name: "ABL1", g2chr: 22, g2start: 23523148, g2end: 23523148, g2name: "BCR"}
                                
                                var circos_events = [];
                                circos_events.push(event);
                                
                                for(var i=0; i<response.data.rows.items.length; i++)
                                {
                                        var item = response.data.rows.items[i];
                                        
                                        // Prototype: {fusion: "FGFR3--TACC3", g1chr: 4, g1start: 1795662, g1end: 1808986, g1name: "FGFR3", g2chr: 4, g2start: 1723217, g2end: 1746905, g2name: "TACC3"},
                                        
                                        if(i<10) console.log(item);
                                        
                                        item.fusion = item.g1 + "--" + item.g2;
                                        item.g1name = item.g1;
                                        item.g2name = item.g2;
                                        delete item.g1
                                        delete item.g2
                                        
                                        item.g1start = parseInt(item.g1start);
                                        item.g1end = item.g1start
                                        item.g2start = parseInt(item.g2start);
                                        item.g2end = item.g2start
                                        item.g1chr = item.g1chr.replace("chr", "");
                                        item.g2chr = item.g2chr.replace("chr", "");
                                        
                                        circos_events.push(item);

//                                         circos_events.push({fusion: gene1+"--"+gene2 , g1chr: chrm1, g1start: fusion_point1, g1end: fusion_point1, g1name: gene1, g2chr: chrm2, g2start: fusion_point2, g2end: fusion_point2, g2name: gene2});
                                }                              
                                
                                console.log("CIRCOS EVENTS", circos_events);
                                circos_events = circos_events.slice(0, 30);
                                
                                self.circos.events[2] = circos_events;
                                
                                BioCircos01 = new BioCircos(self.circos.events, self.circos.genome, {
                                    target : "biocircos",
                                    svgWidth : 600,
                                    svgHeight : 400,
                                    innerRadius: 160,
                                    outerRadius: 180,
                                    genomeFillColor: ["#FFFFCC", "#CCFFFF", "#FFCCCC", "#CCCC99","#0099CC", "#996699", "#336699", "#FFCC33","#66CC00"],
                                    LINKMouseEvent : true,
                                    LINKMouseClickDisplay : true,
                                    LINKMouseClickOpacity : 1.0,
                                    LINKMouseClickStrokeColor : "red",
                                    LINKMouseClickStrokeWidth : 6,
                                    LINKLabelDragEvent : false,
                                    });
                                    BioCircos01.draw_genome(BioCircos01.genomeLength);
                                    
                                }, function(response){
                                        self.sending = false;
                                        console.log("ERROR WHILE GETTING CIRCOS DATA...", response);
                                });
                };
                
                self.load_circos = function(url){
                        $http.get(url).then(function(response) {
                                console.log("CIRCOS DATA OPTIONS: SUCCESS", response.data);
                                self.circos.values = response.data;
                        }, function(response){
                                console.log("ERROR WHILE GETTING CIRCOS DATA OPTIONS...", response);
                        });
                };
                
                self.send_query = function(){
                        var args = [];
                        console.log(self.form)
                        for(var i=0; i<self.form.fields.length; i++)
                        {
                                var field = self.form.fields[i];
                                var value = field.value;
                                console.log("VALUE", value)
                                // if (i == 0 && (value == undefined || value == "")) continue;
                                if (value == undefined || value == "undefined" || value == "") value = "ALL";
                                args.push(value);
                        }
                        
                        console.log("ARGS", args);
                        
                        self.sending = true;
                        if( ! self.form.submit.url.endsWith("/") ) self.form.submit.url = self.form.submit.url  + "/";
                        $http.get(self.form.submit.url + args.join("/") + "/").then(function(response) {
                                self.sending = false;
                                console.log("QUERY SUCCESS", response);
                                self.form_results = response.data.rows;

                                // console.log("HEADER", self.form_results.header);
                                
                        }, function(response){
                                self.sending = false;
                                console.log("ERROR WHILE SENDING QUERY...", response);
                        });
                };
                
                self.ajax2forms = [];

                self.load_data = function(url){
                        
                        $scope.data = [];
                        self.sending = true;
                        $http.get(url).then(function(response)
                {
                                                self.sending = false;
                                                console.log("SUCCESS IN GETTING DATA FROM " + url);
                                                console.log(response.data);
                                                $scope.data = response.data.details;
                                        },
                                        function myError(response) {
                                                self.sending = false;
                                                console.log(response);
                                console.log("ERROR IN GETTING DATA FROM " + url);
                                        }
                        );
                };
                
                self.load_subdata = function(url){
                        
                        $scope.subdata = [];
                        self.sending = true;
                        $http.get(url).then(function(response)
                {
                                                self.sending = false;
                                                console.log("SUCCESS IN GETTING DATA FROM " + url);
                                                console.log(response.data);
                                                $scope.subdata = response.data.details;
                                        },
                                        function myError(response) {
                                                self.sending = false;
                                                console.log(response);
                                console.log("ERROR IN GETTING DATA FROM " + url);
                                        }
                        );
                };
                
                $http.get('config.json').then(function(response) {
                    self.info = response.data;
                    self.info.image.percentage_width_original = self.info.image.percentage_width;
                    console.log(self.info);

                    self.load_data(self.info.links.statistics_all);
                    self.clicked_chromosome = "1"
                    self.load_subdata(self.info.links.statistics_single_chromosome + self.clicked_chromosome + "/");

                    $http.get(self.info.links.statistics_by_chromosome).then(function(response)
                                    {
                                                                    console.log("SUCCESS IN GETTING STATISTICS BY CHROMOSOME");
                                                                    console.log(response.data);
                                                                    $scope.pielabels = response.data.details.header; 
                                                                    $scope.piedata = response.data.details.items;
                                                                    $scope.pieheader = response.data.details.labels;
                                                            },
                                                            function myError(response) {
                                                                    console.log(response);
                                                    console.log("ERROR IN GETTING STATISTICS CHROMOSOME");
                                                            }
                                            );

                    for(var f=0; f<self.info.forms.length; f++)
                    {
                            var form = self.info.forms[f];
                            
                            for(var i=0; i<form.fields.length; i++)
                            {
                                    var field = form.fields[i];
                                    
                                    // Assign the chosen default value
                                    if(field.default) field.value = field.default;
                                    
                                    // Check types
                                    if(field.type == "number") field.value = parseInt(field.value);
                                    
                                    if(field.type === "checkbox" && field.value === "true")
                                    {
                                            field.checked = true;
                                    }
                                    else if(field.type === "select" || (field.type === "checkbox" && field.values.startsWith("http")))
                                    {
                                            var ajaxForms = self.ajax2forms[field.values];
                                            if( !ajaxForms ) {ajaxForms = []; self.ajax2forms[field.values] = ajaxForms;}
                                            ajaxForms.push(field);
                                    }
                                    
                                    // Subform handling 
                                    if(field.form) {
                                            // console.log("Handling subforms");
                                            for(var j=0; j < field.form.length; j++)
                                            {	
                                                    var subform = field.form[j];
                    // 	                        			console.log("SUBFORM: ", subform);
                                                    
                                                    for(var s=0; s<subform.fields.length; s++)
                                                    {
                                                            var subfield = subform.fields[s];
                    // 	                        				console.log("\tSUBFORM FIELD", subfield);
                                                            if(subfield.type == "select") {
                                                                    console.log("FOUND SUBSELECT", subfield);
                                                                    var ajaxForms = self.ajax2forms[subfield.values];
                                                                    if( !ajaxForms ) {ajaxForms = []; self.ajax2forms[subfield.values] = ajaxForms;}
                                                                    ajaxForms.push(subfield);
                                                            }
                                                    }
                                                    
                                                    // console.log("Handling subform = ", subform);
                                                    
                                                    if(!subform.value && !subform.index) continue;
                                                    
                                                    // console.log("\tsubform= ", subform.value, subform.index, field.values);
                                                    
                                                    for(var c=0; c < field.values.length; c++)
                                                            {
                                                                    var value = field.values[c];
                                                                    // console.log("\t\tvalue=", value);
                                                                    
                                                                    if(subform.value && value.label == subform.value)
                                                                    {
                                                                            // console.log("\tsubform's name="+subform.value+" equals name of option="+value);
                                                                            value.url = subform.name;
                                                                            value.form = subform;
                                                                    }
                                                                    else if(subform.index && c == subform.index)
                                                                    {
                                                                            value.url = subform.name;
                                                                            value.form = subform;
                                                                    }
                                                            }
                                            }
                                    }
                            }
                    }

                    console.log("ajax2forms", self.ajax2forms);

                    // Make all the ajax calls here
                    for(var key in self.ajax2forms)
                    {
                            console.log("GETTING VALUES FOR SELECT: " + key);	                        		
                            $http.get(key).then(
                                            function(response) {
                                                    console.log("RESPONSE", response);
                                                    console.log("FORMS AFFECTED", self.ajax2forms[response.config.url]);
                    //    	                        			console.log("ajax2forms", self.ajax2forms);
                    //    	                        			console.log("self", self);
                    //    	                        			console.log("url", response.config.url);
                    //    	                        			console.log("form chosen", self.ajax2forms[response.config.url]);

                                                                            var affectedForms = self.ajax2forms[response.config.url]
                                                                            for(var formId in affectedForms)
                                                                                    affectedForms[formId].values = response.data;
                    //    	                        			console.log("form chosen after", self.ajax2forms[response.config.url]);
                                            
                                            // var selectField = self.form.fields[self.formIndex];	                        			
                                            // selectField.values = response.data;
                                            // console.log("NEW VALUES A: ", selectField);
                                            // console.log("NEW VALUES B: ", self.info.forms.fields[selectIndex]);
                                    },
                                    function(response){
                                            console.log("COULD NOT GET VALUES FOR SELECT FIELD", response);
                                            var affectedForms = self.ajax2forms[response.config.url]
                                                                            for(var formId in affectedForms)
                                                                                    affectedForms[formId].values = ["NOT AVAILABLE"];
                                    }
                        );
                    }

                    self.show_circos(self.info.links.circos_url + self.circos.value + "/");

                    }, function(response){
                            console.log("ERROR GETTING THE CONFIG FILE");
                    });
                
                    self.show_dialog = function(ev, card){
                        
                            console.log(ev, card);
                        
                            $mdDialog.show({
                                locals: {initial_data: card},
                                controller: 'cardController',
                                templateUrl: 'components/card/card.html',
                                parent: angular.element(document.body),
                                targetEvent: ev,
                                clickOutsideToClose:true
                            });
                        };
                                    
                            
        
                    self.exists = function(item, field){
                        if(field.value == undefined) return false;
                        return field.value.indexOf(item) > -1;
                    };
                    self.toggle = function (item, field) {
                        if(field.value == undefined) field.value = []                        
                        var idx = field.value.indexOf(item);
                        if (idx > -1) field.value.splice(idx, 1);
                        else field.value.push(item);
                    };
                                    
                                $scope.chartlabels = ["January", "February", "March", "April", "May", "June", "July"];
                                $scope.chartseries = ['Series A', 'Series B'];
                                $scope.chartdata = [
                                    [65, 59, 80, 81, 56, 55, 40],
                                    [28, 48, 40, 19, 86, 27, 90]
                                ];
                                $scope.onClick = function (points, evt) {
                                    console.log(points, evt);
                                };
                                $scope.datasetOverride = [{ yAxisID: 'y-axis-1' }, { yAxisID: 'y-axis-2' }];
                                $scope.options = {
                                    scales: {
                                        yAxes: [
                                            {
                                                id: 'y-axis-1',
                                                type: 'linear',
                                                display: true,
                                                position: 'left'
                                            },
                                            {
                                                id: 'y-axis-2',
                                                type: 'linear',
                                                display: true,
                                                position: 'right'
                                            }
                                        ]
                                    }
                                };

                                $scope.pieclick = function(evt)
                                {
                                        var chr = $scope.pielabels[evt[0]._index];
                                        var chr_for_query = chr.replace("chr", "").trim();
                                        self.clicked_chromosome = chr_for_query;
                                        self.load_subdata(self.info.links.statistics_single_chromosome + chr_for_query + "/");
                                };
                                $scope.pieoptions = {fontSize: 200};


                                this.settings = {
                                        printLayout: true,
                                        showRuler: true,
                                        showSpellingSuggestions: true,
                                        presentationMode: 'edit'
                                    };

                                    this.sampleAction = function(name, ev) {
                                        $mdDialog.show($mdDialog.alert()
                                            .title(name)
                                            .textContent('You triggered the "' + name + '" action')
                                            .ok('Great')
                                            .targetEvent(ev)
                                        );
                                    };
                            }// End of controller
);