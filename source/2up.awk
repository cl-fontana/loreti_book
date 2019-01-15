/\{geometry\}/ { sub("centering", "hmarginratio=3:5"); }
               { print; }
/frontmatter/  { print "\\null\\clearpage"; }
