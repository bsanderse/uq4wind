function Y = NM80_readoutput(output_dir,P)


% return output depending on the QoI we are interested in
switch P{27}
    case 'Power'
        filename = fullfile(output_dir,'AeroPower.dat'); % location of the AeroPower.dat output file
        [Times,Azimuthdeg,PowerWatt,Axial_ForceN] = AeroPower(filename, P{23}, P{24});
        Y = mean(PowerWatt);
    case 'Axial_Force'
        filename = fullfile(output_dir,'AeroPower.dat'); % location of the AeroPower.dat output file
        [Times,Azimuthdeg,PowerWatt,Axial_ForceN] = AeroPower(filename, P{23}, P{24});
        Y = mean(Axial_ForceN);
        
    otherwise
        error(strcat('QoI type unknown; check the turbine file ',P{29}));
end

end