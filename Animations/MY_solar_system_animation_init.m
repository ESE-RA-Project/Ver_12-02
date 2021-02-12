% time_vector, a row for each day of the mission
time_vector = ymd_gen([2022, 8, 2],[2028, 5, 29]);
n_days = size(time_vector,1);

%% Mission info
%{
    The mission is subdivided in three parts:
    - Earth to Venus;
    - Venus to Earth;
    - Earth to Saturn;
    - Saturn to Enceladus.
    Variables for each part have the same name but different suffix.
    
    Planet configurations have the following numerical indications:
    - 0: at Earth departure (to Venus) -> 02/08/2022
    - 1: at Venus arrival/departure (for the first fly-by) -> vedere data 
    - 2: at Earth arrival/departure (for the second fly-by) -> per i flyby
    - 3: at Saturn arrival (from Earth) -> 02/03/2030
    - 4: at Saturn departure (to Enceladus) -> vedere quando sono + vicini 
    - 5: at Enceladus arrival (from Saturn) -> saturno ed Enceladus
%}

%% Plots Parameters

% init of status msg
status_msg = ['Nothing'];

% axis lim in plots
xy_lim = 1.5e9;		%lim in xy plane
z_lim = 5e8;		%lim in z coord

% colours of objects, rgb
% nice link to get rbg colour:
% https://www.rapidtables.com/web/color/RGB_Color.html
% col = [	"g"		     %green
% 		"m"          %magenta
% 		"b"          %blue
% 		"r"          %red
%         "#DFEC26"     %giallo venere
% 		"#A2142F"    %darker red
% 		"#7E2F8E"    %purple
% 		"#4DBEEE"    %darker cyan
% 		"c"          %(bright) cyan
% 		"#ECF38D"    %giallo saturno
% 		"#70FF9B"    %darker green
% 		"#EDB120"    %ochre
% 		"#D95319"];  %orange, not visible due to Sun orbit dimensions
col = ["g"          %green
      "m"          %magenta
      "b"          %blue
      "r"          %red
      "#A2142F"    %darker red
      "#7E2F8E"    %purple
      "#4DBEEE"    %darker cyan
      "c"          %(bright) cyan
      "#D95319"    %orange
      "#77AC30"    %darker green
      "#D95319"];  %orange, not visible due to Sun orbit dimensions

col_bkgnd	= [0,	0,	0];		% Background
col_grid	= [255,	255, 255]	/255;	% Grid
col_sun		= [255, 204, 0]		/255;	
% col_sun=col(11);
col_earth	= col(3);	%[0,	102, 204]	/255;
col_venus	= col(2);	%[255,	102, 0]		/255;	
col_saturn	= col(6);	%[0,	102, 204]	/255;	
col_enceladus	= col(10);	%[102,	153, 153]	/255;	
col_spcr	= [0,	196, 255]	/255;	% spacecraft
col_jupiter=col(5);
col_mars=col(4);

% dimension of objects
dim_sun		= 20;
dim_earth	= 8;
dim_venus	= 7;
dim_saturn	= 15;
dim_jupiter	= 17;
dim_enceladus	= 5;
dim_spcr	= 3;
dim_mars = 6;

% width of lines
planet_linewidth = 0.1;
width_spcr  = 1;

%% Movie parameters

% init of frame counter
k = 1;
movie_fps = 20;