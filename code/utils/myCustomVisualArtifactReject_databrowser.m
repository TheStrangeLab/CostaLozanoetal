function [newdata, trialsrejected] = myCustomVisualArtifactReject_databrowser(cfg,data)
samples = data.sampleinfo;
cfg_new = ft_databrowser(cfg,data);
cfg_new.artfctdef.reject = 'complete';
newdata = ft_rejectartifact(cfg_new, data)
artifacts = newdata.cfg.artfctdef.visual.artifact;

z=1;
for s = 1:size(samples,1)
    for a = 1:size(artifacts,1)
        if (artifacts(a,1) >= samples(s,1)) && (artifacts(a,2) <= samples(s,2))
        trialsrejected(z) = s;
        z = z + 1;
        end
    end
end

pause(0.1);

%[TF,trialsrejected] = ismember(artifacts,samples,'rows');% [];
%I have to check what trial correspond to each artifact
end