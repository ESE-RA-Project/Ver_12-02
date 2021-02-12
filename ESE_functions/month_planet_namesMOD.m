function [month, planet] = month_planet_names(month_id, planet_id)
% MONTH_PLANET_NAMES returns the name of the month and the planet
%   corresponding, respectively, to the numbers "month_id" and
%   "planet_id".
%  
%   months    - a vector containing the names of the 12 months
%   planets   - a vector containing the names of the 9 planets, plus
%               Enceladus and the Sun
%   month_id  - the month number (1 - 12)
%   planet_id - the planet number (1 - 11)
%  
% User M-functions required: none

    months  = ['January  '
               'February '
               'March    '
               'April    '
               'May      '
               'June     '
               'July     '
               'August   '
               'September'
               'October  '
               'November '
               'December '];

    planets = ['Mercury  '
               'Venus    '
               'Earth    '
               'Mars     '
               'Jupiter  '
               'Saturn   '
               'Uranus   '
               'Neptune  '
               'Pluto    '
               'Enceladus'
               'Sun      '];

    month   = months (month_id,  1:9);
    planet  = planets(planet_id, 1:9);

end %month_planet_names
