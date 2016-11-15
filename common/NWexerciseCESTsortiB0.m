function out = NWexerciseCESTsortiB0(in)
%
%
% Guesses scan directories in exercise CEST protocol
%
% Use 'pre' and 'post' for pre and post exercise
% Use 'ref' or 'none' for reference scan


if nargin>0
    scandirs = dir(in);
else
    scandirs = dir;
end

out = []; 
for zz = 1:2 % loop through 2 times
    for ii=1:length(scandirs)
        if scandirs(ii).isdir
            lowname = lower(scandirs(ii).name);
            splitname1 = strsplit(lowname,'_');
            lowname1 = strjoin(splitname1(3:end),'_'); % ignores sequence name part
            if strfind(lowname,'b0') % B0
                if strfind(lowname,'post')
                    splitname = strsplit(lowname,'_');
                    if strcmp(lowname(1),'s') && ~isnan(str2double(lowname(2:5)))
                        ind = str2double(lowname(2:5));
                    else
                        ind = str2double(splitname(end));
                    end
                    if zz==1
                        if exist('B0post1','var')
                            B0post2 = ind;
                        else
                            B0post1 = ind;
                        end
                    elseif isequal(B0post1+1,B0post2)
                        if isequal(ind,B0post1)
                            out.Dicom.B0magpost = scandirs(ii).name;
                        elseif isequal(ind,B0post2)
                            out.Dicom.B0phpost = scandirs(ii).name;
                        else
                            disp([ii,zz])
                            warning('something unexpected')
                        end
                    elseif isequal(B0post2+1,B0post1)
                        if isequal(ind,B0post2)
                            out.Dicom.B0magpost = scandirs(ii).name;
                        elseif isequal(ind,B0post1)
                            out.Dicom.B0phpost = scandirs(ii).name;
                        else
                            disp([ii,zz])
                            warning('something unexpected')
                        end
                    else
                        disp([ii,zz])
                        warning('something unexpected')
                    end
                elseif strfind(lowname,'_pre')
                    splitname = strsplit(lowname,'_');
                    if strcmp(lowname(1),'s') && ~isnan(str2double(lowname(2:5)))
                        ind = str2double(lowname(2:5));
                    else
                        ind = str2double(splitname(end));
                    end
                    if zz==1
                        if exist('B0pre1','var')
                            B0pre2 = ind;
                        else
                            B0pre1 = ind;
                        end
                    elseif isequal(B0pre1+1,B0pre2)
                        if isequal(ind,B0pre1)
                            out.Dicom.B0magpre = scandirs(ii).name;
                        elseif isequal(ind,B0pre2)
                            out.Dicom.B0phpre = scandirs(ii).name;
                        else
                            disp([ii,zz])
                            warning('something unexpected')
                        end
                    elseif isequal(B0pre2+1,B0pre1)
                        if isequal(ind,B0pre2)
                            out.Dicom.B0magpre = scandirs(ii).name;
                        elseif isequal(ind,B0pre1)
                            out.Dicom.B0phpre = scandirs(ii).name;
                        else
                            disp([ii,zz])
                            warning('something unexpected')
                        end
                    else
                        disp([ii,zz])
                        warning('something unexpected')
                    end
                else
                    splitname = strsplit(lowname,'_');
                    if strcmp(lowname(1),'s') && ~isnan(str2double(lowname(2:5)))
                        ind = str2double(lowname(2:5));
                    else
                        ind = str2double(splitname(end));
                    end
                    if zz==1
                        if exist('B0unk1','var')
                            B0unk2 = ind;
                        else
                            B0unk1 = ind;
                        end
                    else
                        if isequal(B0unk1+1,B0unk2)
                            if ~isfield(out.Dicom,'B0magpost') && isfield(out.Dicom,'B0magpre')
                                if isequal(ind,B0unk1)
                                    out.Dicom.B0magpost = scandirs(ii).name;
                                elseif isequal(ind,B0unk2)
                                    out.Dicom.B0phpost = scandirs(ii).name;
                                else
                                    disp([ii,zz])
                                    warning('something unexpected')
                                end
                            elseif ~isfield(out.Dicom,'B0pre') && isfield(out.Dicom,'B0post')
                                if isequal(ind,B0unk1)
                                    out.Dicom.B0magpre = scandirs(ii).name;
                                elseif isequal(ind,B0unk2)
                                    out.Dicom.B0phpre = scandirs(ii).name;
                                else
                                    disp([ii,zz])
                                    warning('something unexpected')
                                end
                            end
                        elseif isequal(B0unk2+1,B0unk1)
                            if ~isfield(out.Dicom,'B0magpost') && isfield(out.Dicom,'B0magpre')
                                if isequal(ind,B0unk2)
                                    out.Dicom.B0magpost = scandirs(ii).name;
                                elseif isequal(ind,B0unk1)
                                    out.Dicom.B0phpost = scandirs(ii).name;
                                else
                                    disp([ii,zz])
                                    warning('something unexpected')
                                end
                            elseif ~isfield(out.Dicom,'B0pre') && isfield(out.Dicom,'B0post')
                                if isequal(ind,B0unk2)
                                    out.Dicom.B0magpre = scandirs(ii).name;
                                elseif isequal(ind,B0unk1)
                                    out.Dicom.B0phpre = scandirs(ii).name;
                                else
                                    disp([ii,zz])
                                    warning('something unexpected')
                                end
                            end
                        else
                            disp([ii,zz])
                            warning('something unexpected')
                        end
                    end
                end
            elseif strfind(lowname1,'b1') % B1
                if strfind(lowname,'post')
                    out.Dicom.B1post = scandirs(ii).name;
                elseif strfind(lowname,'_pre')
                    out.Dicom.B1pre = scandirs(ii).name;
                else
                    if zz==2
                        if ~isfield(out.Dicom,'B1post') && isfield(out.Dicom,'B1pre')
                            out.Dicom.B1post = scandirs(ii).name;
                        elseif ~isfield(out.Dicom,'B1pre') && isfield(out.Dicom,'B1post')
                            out.Dicom.B1pre = scandirs(ii).name;
                        else
                            out.Dicom.unknown.B1 = scandirs(ii).name;
                        end
                    end
                end
            elseif strfind(lowname1,'cest') % CEST images
                if strfind(lowname,'post')
                    out.Dicom.CESTpost = scandirs(ii).name;
                elseif strfind(lowname,'_pre')
                    out.Dicom.CESTpre = scandirs(ii).name;
                else
                    if zz==2
                        if ~isfield(out.Dicom,'CESTpost') && isfield(out.Dicom,'CESTpre')
                            out.Dicom.CESTpost = scandirs(ii).name;
                        elseif ~isfield(out.Dicom,'CESTpre') && isfield(out.Dicom,'CESTpost')
                            out.Dicom.CESTpre = scandirs(ii).name;
                        else
                            out.Dicom.unknown.CEST = scandirs(ii).name;
                        end
                    end
                end
%                 if strfind(lowname,'post')
%                     splitname = strsplit(lowname,'_');
%                     if strcmp(lowname(1),'s') && ~isnan(str2double(lowname(2:5)))
%                         ind = str2double(lowname(2:5));
%                     else
%                         ind = str2double(splitname(end));
%                     end
%                     if zz==1
%                         if exist('CESTpost1','var')
%                             CESTpost2 = ind;
%                         else
%                             CESTpost1 = ind;
%                         end
%                     elseif isequal(CESTpost1+1,CESTpost2)
%                         if isequal(ind,CESTpost1)
%                             out.Dicom.CESTpost = scandirs(ii).name;
%                         end
%                     elseif isequal(CESTpost2+1,CESTpost1)
%                         if isequal(ind,CESTpost2)
%                             out.Dicom.CESTpost = scandirs(ii).name;
%                         end
%                     else
%                         disp([ii,zz])
%                         warning('something unexpected')
%                     end
%                 elseif strfind(lowname,'_pre')
%                     splitname = strsplit(lowname,'_');
%                     if strcmp(lowname(1),'s') && ~isnan(str2double(lowname(2:5)))
%                         ind = str2double(lowname(2:5));
%                     else
%                         ind = str2double(splitname(end));
%                     end
%                     if zz==1
%                         if exist('CESTpre1','var')
%                             CESTpre2 = ind;
%                         else
%                             CESTpre1 = ind;
%                         end
%                     elseif isequal(CESTpre1+1,CESTpre2)
%                         if isequal(ind,CESTpre1)
%                             out.Dicom.CESTpre = scandirs(ii).name;
%                         end
%                     elseif isequal(CESTpre2+1,CESTpre1)
%                         if isequal(ind,CESTpre2)
%                             out.Dicom.CESTpre = scandirs(ii).name;
%                         end
%                     else
%                         disp([ii,zz])
%                         warning('something unexpected')
%                     end
%                 else
%                     splitname = strsplit(lowname,'_');
%                     if strcmp(lowname(1),'s') && ~isnan(str2double(lowname(2:5)))
%                         ind = str2double(lowname(2:5));
%                     else
%                         ind = str2double(splitname(end));
%                     end
%                     if zz==1
%                         if exist('B0unk1','var')
%                             CESTunk2 = ind;
%                         else
%                             CESTunk1 = ind;
%                         end
%                     else
%                         if isequal(CESTunk1+1,CESTunk2)
%                             if ~isfield(out.Dicom,'CESTpost') && isfield(out.Dicom,'CESTpre')
%                                 if isequal(ind,CESTunk1)
%                                     out.Dicom.CESTpost = scandirs(ii).name;
%                                 else
%                                     disp([ii,zz])
%                                     warning('something unexpected')
%                                 end
%                             elseif ~isfield(out.Dicom,'CESTpre') && isfield(out.Dicom,'CESTpost')
%                                 if isequal(ind,CESTunk1)
%                                     out.Dicom.CESTpre = scandirs(ii).name;
%                                 else
%                                     disp([ii,zz])
%                                     warning('something unexpected')
%                                 end
%                             end
%                         elseif isequal(CESTunk2+1,CESTunk1)
%                             if ~isfield(out.Dicom,'CESTpost') && isfield(out.Dicom,'CESTpre')
%                                 if isequal(ind,CESTunk2)
%                                     out.Dicom.CESTpost = scandirs(ii).name;
%                                 else
%                                     disp([ii,zz])
%                                     warning('something unexpected')
%                                 end
%                             elseif ~isfield(out.Dicom,'CESTpre') && isfield(out.Dicom,'CESTpost')
%                                 if isequal(ind,CESTunk2)
%                                     out.Dicom.CESTpre = scandirs(ii).name;
%                                 else
%                                     disp([ii,zz])
%                                     warning('something unexpected')
%                                 end
%                             end
%                         else
%                             disp([ii,zz])
%                             warning('something unexpected')
%                         end
%                     end
%                 end
            elseif strfind(lowname1,'wassr') % WASSR images
                if strfind(lowname,'post')
                    out.Dicom.WASSRpost = scandirs(ii).name;
                elseif strfind(lowname,'_pre')
                    out.Dicom.WASSRpre = scandirs(ii).name;
                else
                    if zz==2
                        if ~isfield(out.Dicom,'WASSRpost') && isfield(out.Dicom,'WASSRpre')
                            out.Dicom.WASSRpost = scandirs(ii).name;
                        elseif ~isfield(out.Dicom,'WASSRpre') && isfield(out.Dicom,'WASSRpost')
                            out.Dicom.WASSRpre = scandirs(ii).name;
                        else
                            out.Dicom.unknown.WASSR = scandirs(ii).name;
                        end
                    end
                end
            elseif ~(isempty(strfind(lowname1,'none')) && isempty(strfind(lowname1,'ref'))) % Reference scan
                if strfind(lowname,'post')
                    out.Dicom.Refpost = scandirs(ii).name;
                elseif strfind(lowname,'_pre')
                    out.Dicom.Refpre = scandirs(ii).name;
                else
                    if zz==2
                        if ~isfield(out.Dicom,'Refpost') && isfield(out.Dicom,'Refpre')
                            out.Dicom.Refpost = scandirs(ii).name;
                        elseif ~isfield(out.Dicom,'Refpre') && isfield(out.Dicom,'Refpost')
                            out.Dicom.Refpre = scandirs(ii).name;
                        else
                            out.Dicom.unknown.Ref = scandirs(ii).name;
                        end
                    end
                end
            end
        elseif strcmp(scandirs(ii).name(end-3:end),'.dat')
            splitname1 = strsplit(scandirs(ii).name,'_');
            lowname = lower(strjoin(splitname1(6:end),'_')); % ignores sequence name part
            out.Raw.unknown{1} = [];
            if ~(isempty(strfind(lowname,'none'))) || ~(isempty(strfind(lowname,'ref')))
                if ~(isempty(strfind(lowname,'pre')))
                    out.Raw.Refpre = scandirs(ii).name;
                elseif ~isempty(strfind(lowname,'post'))
                    out.Raw.Refpost = scandirs(ii).name;
                else
                    out.Raw.unknown{end+1} = scandirs(ii).name;
                end
            elseif ~(isempty(strfind(lowname,'b1')))
                if ~(isempty(strfind(lowname,'pre')))
                    out.Raw.B1pre = scandirs(ii).name;
                elseif ~isempty(strfind(lowname,'post'))
                    out.Raw.B1post = scandirs(ii).name;
                else
                    out.Raw.unknown{end+1} = scandirs(ii).name;
                end
            elseif ~(isempty(strfind(lowname,'b0')))
                if ~(isempty(strfind(lowname,'pre')))
                    out.Raw.B0pre = scandirs(ii).name;
                elseif ~isempty(strfind(lowname,'post'))
                    out.Raw.B0post = scandirs(ii).name;
                else
                    out.Raw.unknown{end+1} = scandirs(ii).name;
                end
            elseif ~(isempty(strfind(lowname,'cest')))
                if ~(isempty(strfind(lowname,'pre')))
                    out.Raw.CESTpre = scandirs(ii).name;
                elseif ~isempty(strfind(lowname,'post'))
                    out.Raw.CESTpost = scandirs(ii).name;
                else
                    out.Raw.unknown{end+1} = scandirs(ii).name;
                end
            elseif ~(isempty(strfind(lowname,'iwassr')))
                if ~(isempty(strfind(lowname,'pre')))
                    out.Raw.iWASSRpre = scandirs(ii).name;
                elseif ~isempty(strfind(lowname,'post'))
                    out.Raw.iWASSRpost = scandirs(ii).name;
                else
                    out.Raw.unknown{end+1} = scandirs(ii).name;
                end
            elseif ~(isempty(strfind(lowname,'wassr')))
                if ~(isempty(strfind(lowname,'pre')))
                    out.Raw.WASSRpre = scandirs(ii).name;
                elseif ~isempty(strfind(lowname,'post'))
                    out.Raw.WASSRpost = scandirs(ii).name;
                else
                    out.Raw.unknown{end+1} = scandirs(ii).name;
                end
            end
        end
    end
end

if isfield(out,'Dicom')
    if ~isfield(out.Dicom,'Refpost'), out.Dicom.Refpost = ''; end
    if ~isfield(out.Dicom,'Refpre'), out.Dicom.Refpre = ''; end
    if ~isfield(out.Dicom,'WASSRpost'), out.Dicom.WASSRpost = ''; end
    if ~isfield(out.Dicom,'WASSRpre'), out.Dicom.WASSRpre = ''; end
    if ~isfield(out.Dicom,'CESTpost'), out.Dicom.CESTpost = ''; end
    if ~isfield(out.Dicom,'CESTpre'), out.Dicom.CESTpre = ''; end
    if ~isfield(out.Dicom,'B1post'), out.Dicom.B1post = ''; end
    if ~isfield(out.Dicom,'B1pre'), out.Dicom.B1pre = ''; end
    if ~isfield(out.Dicom,'Refpost'), out.Dicom.Refpost = ''; end
    if ~isfield(out.Dicom,'B0magpost') || ~isfield(out.Dicom,'B0phpost'), out.Dicom.B0magpost = ''; out.Dicom.B0phpost = ''; end
    if ~isfield(out.Dicom,'B0magpre') || ~isfield(out.Dicom,'B0phpre'), out.Dicom.B0magpre = ''; out.Dicom.B0phpre = ''; end
end

if isfield(out,'Raw')
    if ~isfield(out.Raw,'Refpre'),out.Raw.Refpre = ''; end
    if ~isfield(out.Raw,'Refpost'),out.Raw.Refpost = ''; end
    if ~isfield(out.Raw,'B1pre'),out.Raw.B1pre = ''; end
    if ~isfield(out.Raw,'B1post'),out.Raw.B1post = ''; end
    if ~isfield(out.Raw,'B0pre'),out.Raw.B0pre = ''; end
    if ~isfield(out.Raw,'B0post'),out.Raw.B0post = ''; end
    if ~isfield(out.Raw,'CESTpre'),out.Raw.CESTpre = ''; end
    if ~isfield(out.Raw,'CESTpost'),out.Raw.CESTpost = ''; end
    if ~isfield(out.Raw,'WASSRpre'),out.Raw.WASSRpre = ''; end
    if ~isfield(out.Raw,'WASSRpost'),out.Raw.WASSRpost = ''; end
end

