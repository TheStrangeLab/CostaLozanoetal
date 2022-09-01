function [eog,perif] = getcapinfo(captype)

	if ~exist('captype','var')
		captype = 'GSN200';
	end
	
	switch captype
		case 'GSN200'
		eog = {[26 127], [8 126]};
		perif = [127 126 17 128 125 120 44 49 56 63 69 74 82 89 95 100 108 114];
		case 'HCGSN'
		eog = {[25 127], [8 126]};
		perif = [127 126 17 21 14 25 8 128 125 48 119 49 113];
	end
