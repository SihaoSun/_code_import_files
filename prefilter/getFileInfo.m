

function [FileNames, FileCount, Flag] = getFileInfo(fnames, path)

    FileCount = 0;
    FileNames = [];
    Flag = 1;
    
    if (iscell(fnames))
        for j = 1:length(fnames)
            found = false;
            fname = strcat(path, fnames{j});
            for i = 1:size(FileNames, 1)
                if (strcmpi(fname, FileNames{i, 1}))
                    found = true;
                    break;
                end
            end
            if (~found)
                try
                    finfo = dir(fname);
                catch
                    warning('Could not load data from <%s>, returning...', fname);
                    Flag = 0;
                    return;
                end
                if (isempty(finfo))
                    warning('Could not load data from <%s>, returning...', fname);
                    Flag = 0;
                    return;
                end
                FileCount = FileCount + 1;
                FileNames{FileCount, 1} = fname;
                FileNames{FileCount, 2} = fnames{j};
                FileNames{FileCount, 3} = finfo.bytes;
            end
        end
    else
        found = false;
        fname = strcat(path, fnames);
        for i = 1:size(FileNames, 1)
            if (strcmpi(fname, FileNames{i, 1}))
                found = true;
                break;
            end
        end
        if (~found)
            try
                finfo = dir(fname);
            catch
                warning('Could not load data from <%s>, returning...', fname);
                Flag = 0;
                return;
            end
            if (isempty(finfo))
                warning('Could not load data from <%s>, returning...', fname);
                Flag = 0;
                return;
            end
            FileCount = FileCount + 1;
            FileNames{FileCount, 1} = fname;
            FileNames{FileCount, 2} = fnames;
            FileNames{FileCount, 3} = finfo.bytes;
        end
    end


