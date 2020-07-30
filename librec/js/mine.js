$(document).ready(function(){    
	// jQuery for page scrolling feature - requires jQuery Easing plugin
	$('a.page-scroll').bind('click', function(event) {
		var entry = $(this);
		$('html, body').stop().animate({
			scrollTop: ($(entry.attr('href')).offset().top-55)
		}, 1000);
		event.preventDefault();
	});	
});