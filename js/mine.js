jQuery.fn.dataTableExt.oSort['shortdate-asc']  = function(x,y) {
	var xi=x.indexOf("<br>");
	if(xi>-1) x = x.substr(0, xi)
	xi=x.indexOf("-");
	if(xi>-1) {
		var x1 = x.substr(0, xi);
		var x2 = x.substr(x.indexOf(","));
		x=x1+x2;
	}
	
	
	var yi=y.indexOf("<br>");
	if(yi>-1) y = y.substr(0, yi)
	yi=y.indexOf("-");
	if(yi>-1) {
		var y1 = y.substr(0, yi);
		var y2 = y.substr(y.indexOf(","));
		y=y1+y2;								
	}
	
	x=Date.parse(x);
	y=Date.parse(y);
	return ((x < y) ? -1 : ((x > y) ?  1 : 0));
};


jQuery.fn.dataTableExt.oSort['shortdate-desc']  = function(x,y) {
	var xi=x.indexOf("<br>");
	if(xi>-1) x = x.substr(0, xi)
	xi=x.indexOf("-");
	if(xi>-1) {
		var x1 = x.substr(0, xi);
		var x2 = x.substr(x.indexOf(","));
		x=x1+x2;
	}
	
	var yi=y.indexOf("<br>");
	if(yi>-1) y = y.substr(0, yi)
	yi=y.indexOf("-");
	if(yi>-1) {
		var y1 = y.substr(0, yi);
		var y2 = y.substr(y.indexOf(","));
		y=y1+y2;								
	}
			 
	x=Date.parse(x);
	y=Date.parse(y);
 
	return ((x < y) ?  1 : ((x > y) ? -1 : 0));
};

$(document).ready(function(){        	
	
	$("#confTable").dataTable({				 
	 "sPaginationType": "full_numbers",
		 "aLengthMenu": [/*real values*/[15, 25, 50, 100, -1], /*display values*/[15, 25, 50, 100, "All"]],
		 'iDisplayLength': 15,
		 "aaSorting": [[3, 'desc'],[1, 'asc']],
		 "bSort":true,
		 "bSortClasses": false,
		 "aoColumnDefs": [
		  { "sSortDataType": "shortdate", "aTargets": [ 2, 3, 4, 5 ] },
		  { "sType": "shortdate", "aTargets": [ 2, 3, 4, 5 ] },
		]
	});
	
	$("#pubTable").dataTable({				 
	 "sPaginationType": "full_numbers",
		 "aLengthMenu": [/*real values*/[10, 25, 50, 100, -1], /*display values*/[10, 25, 50, 100, "All"]],
		 'iDisplayLength': 25,
		 "aoColumns": [null, null, /*{"bVisible": false}, null,*/ null,null],
		 "aaSorting": [[0, 'desc'], [1, 'asc'], [2, 'asc']],
		 "bSort":true,
		 "bSortClasses": false
	});
	
	// jQuery for page scrolling feature - requires jQuery Easing plugin
	$('a.page-scroll').bind('click', function(event) {
		var entry = $(this);
		$('html, body').stop().animate({
			scrollTop: ($(entry.attr('href')).offset().top-50)
		}, 1000);
		event.preventDefault();
	});
});