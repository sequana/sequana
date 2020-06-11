jQuery.tableOfContents =
function (tocList) {
    jQuery(tocList).empty();
    var prevH2Item = null;
    var prevH2List = null;

    var index = 0;
    jQuery("h2, h3").each(function() {

        var anchor = "<a name='" + index + "'></a>";
        jQuery(this).before(anchor);

        var li     = "<li><a href='#" + index + "'>" +
                     jQuery(this).text() + "</a></li>";

        if( jQuery(this).is("h2") ){
            prevH2List = jQuery("<ul></ul>");
            prevH2Item = jQuery(li);
            prevH2Item.append(prevH2List);
            prevH2Item.appendTo(tocList);
        } else {
            prevH2List.append(li);
        }
        index++;
    });
}
