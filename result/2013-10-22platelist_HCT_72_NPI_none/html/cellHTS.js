// switch class labels for tabs
toggleSelectionState = function(id, thisTable)
{
    var tables = document.getElementsByTagName("table");
    var tl = tables.length;
    var cl;
    var match;
    for(var i=0; i<tl; i++)
    {
        cl = tables[[i]].className;
	mt = cl.match(id);
        if(mt)
	    tables[[i]].className = cl.replace(/ selected/, " unselected");
    }
    thisTable.className = thisTable.className.replace(/ unselected/, " selected");
}


// hide or show images of chtsImageStacks
toggleImageVisibility = function(id)
    {
	var tables = document.getElementsByTagName("table");
	var tl = tables.length;
	var cl;
	var mtid; var mtChanRep;
	for(var i=0; i<tl; i++)
	{
            cl = tables[[i]].className;
	    mtId = cl.match(id);
	    mtChanRep = cl.match("channel"+currentChannel+" replicate"+currentReplicate);
	    if(mtId)
		tables[[i]].className = cl.replace(/ visible/, " invisible");
	    if(mtChanRep)
		tables[[i]].className = cl.replace(/ invisible/, " visible");
	}
    }


// Toggle the main table tags
function toggleTabById(id, thisTable, src)
{
    toggleSelectionState(id, thisTable);
    resetIframe("main");
    document.getElementsByTagName("iframe")[0].src=src; 
}



// state variables
var currentReplicate=1;
var currentChannel=1;
var tt_Enabled; // Allows to (temporarily) suppress tooltips, e.g. by providing the user with a button that sets this global variable to false



// Toggle chtsImageStacks by Channel
function toggleTabByChannel(id, thisTable, channel)
{
    toggleSelectionState(id, thisTable);  
    currentChannel = channel;
    id = id.replace(/Channel/, "");
    toggleImageVisibility(id);
}



// Toggle chtsImageStacks by Replicate
function toggleTabByReplicate(id, thisTable, replicate)
{
    toggleSelectionState(id, thisTable);  
    currentReplicate = replicate;
    id = id.replace(/Replicate/, "");
    toggleImageVisibility(id);
}



// Link to pdf version of images
function linkToPdf(url)
{
    document.location.href = url;
    resetIframe("main");
}



// Set current url of the document
function linkToFile(url)
{
    document.location.href = url; 
    parent.resetIframe("main");  
}


// Setter methods for coockies
function Set_Cookie( name, value, expires, path, domain, secure )
{
// set time, it's in milliseconds
var today = new Date();
today.setTime( today.getTime() );

/*
if the expires variable is set, make the correct
expires time, the current script below will set
it for x number of days, to make it for hours,
delete * 24, for minutes, delete * 60 * 24
*/
if ( expires )
{
expires = expires * 1000 * 60 * 60 * 24;
}
var expires_date = new Date( today.getTime() + (expires) );

document.cookie = name + "=" +escape( value ) +
( ( expires ) ? ";expires=" + expires_date.toGMTString() : "" ) +
( ( path ) ? ";path=" + path : "" ) +
( ( domain ) ? ";domain=" + domain : "" ) +
( ( secure ) ? ";secure" : "" );
}



// Getter method for cookies
function Get_Cookie( check_name ) {
	// first we'll split this cookie up into name/value pairs
	// note: document.cookie only returns name=value, not the other components
	var a_all_cookies = document.cookie.split( ';' );
	var a_temp_cookie = '';
	var cookie_name = '';
	var cookie_value = '';
	var b_cookie_found = false; // set boolean t/f default f

	for ( i = 0; i < a_all_cookies.length; i++ )
	{
		// now we'll split apart each name=value pair
		a_temp_cookie = a_all_cookies[i].split( '=' );


		// and trim left/right whitespace while we're at it
		cookie_name = a_temp_cookie[0].replace(/^\s+|\s+$/g, '');

		// if the extracted name matches passed check_name
		if ( cookie_name == check_name )
		{
			b_cookie_found = true;
			// we need to handle case where cookie has no value but exists (no = sign, that is):
			if ( a_temp_cookie.length > 1 )
			{
				cookie_value = unescape( a_temp_cookie[1].replace(/^\s+|\s+$/g, '') );
			}
			// note that in cases where cookie is initialized but no value, null is returned
			return cookie_value;
			break;
		}
		a_temp_cookie = null;
		cookie_name = '';
	}
	if ( !b_cookie_found )
	{
		return null;
	}
}



// Turn tooltips on or off 
function toggleHelp(x)
{
    cn = "showTooltips";
    tt_Enabled = !tt_Enabled
    Set_Cookie(cn, tt_Enabled, '', '/', '', '' );
    if(x.innerHTML == "on") x.innerHTML="off"; else x.innerHTML="on";
}


// Initialize the page
function initialize()
{
    cn = "showTooltips";
    if(!Get_Cookie(cn))
	Set_Cookie(cn, true, '', '/', '', '' );
    tt_Enabled = Get_Cookie(cn);
    cv = Get_Cookie(cn);
    help = document.getElementById("helpSwitch");
    if(cv=="true") 
    {
	tt_Enabled=true; 
	if(help)
	    help.innerHTML="on";
    }
    else 
    {
	tt_Enabled=false;
	if(help)
	    help.innerHTML="off";
    }
   
}


// Automatically resize the iFrame based on its content
function autoIframe(frameId){
    try{
	frame = document.getElementById(frameId);
	innerDoc = (frame.contentDocument) ? frame.contentDocument : frame.contentWindow.document;
	objToResize = (frame.style) ? frame.style : frame;
	objToResize.height = innerDoc.body.scrollHeight + 10;
	objToResize.width = innerDoc.body.scrollWidth + 10;
    }
    catch(err){
	window.status = err.message;
    }
}



// get the current window size (cross browser)
window.size = function(width)
{
	var w = 0;
	var h = 0;

	//IE
	if(!window.innerWidth)
	{
		//strict mode
		if(!(document.documentElement.clientWidth == 0))
		{
			w = document.documentElement.clientWidth;
			h = document.documentElement.clientHeight;
		}
		//quirks mode
		else
		{
			w = document.body.clientWidth;
			h = document.body.clientHeight;
		}
	}
	//w3c
	else
	{
		w = window.innerWidth;
		h = window.innerHeight;
	}
    if(width)
	return w;
    else
	return h
}


// reset the iFrame size to the available window size
function resetIframe(frameId){
    try{
	xoff = 175;
	yoff = 120;
	frame = document.getElementById(frameId);
	objToResize = (frame.style) ? frame.style : frame;
	objToResize.height = window.size(false) - yoff;
	objToResize.width = window.size(true) - xoff;
    }
    catch(err){
	window.status = err.message;
    }
}