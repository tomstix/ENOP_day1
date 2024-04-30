function people(db_file, operation, varargin)
    try 
        db = load(db_file, "db");
        db = db.db;
    catch
        db = table;
        db.Name = {};
        db.Age = {};
    end
    switch operation
        case 'reset'
            % reset database
            db = table;
            db.Name = {};
            db.Age = {};
        case 'list'
            % list all entries
            for i = 1:height(db)
                name = db.Name{i};
                age = db.Age{i};
                fprintf("Name: %s, Age: %d\n", name, age)
            end
        case 'insert'
            % insert new entry
            if nargin > 2 && mod(nargin, 2) ~= 0
                warning("Invalid number of arguments!")
                return
            end
            for i = 1:2:(nargin-2)
                name = varargin{i};
                age = varargin{i+1};
                % Look whether the name already exists
                found_idx = find(db.Name == string(name));
                % variable to determine whether to skip entry
                skip = false;
                if numel(found_idx) > 0
                    for j = 1:numel(found_idx)
                        idx = found_idx(j);
                        if db.Age{idx} == age
                            warning("%s, %d already found!", name, age)
                            skip = true;
                        end
                    end
                end
                if skip
                    continue
                end
                fprintf("Inserting %s, %d\n", name, age)
                cell = {name, age};
                db = [db;cell];
            end
        case 'remove'
            % remove entry
            if nargin > 2 && mod(nargin, 2) ~= 0
                warning("Invalid number of arguments!")
                return
            end
            for i = 1:2:(nargin-2)
                name = varargin{i};
                age = varargin{i+1};
                % Look whether the name already exists
                found_idx = find(db.Name == string(name));
                % variable to determine whether to skip entry
                if numel(found_idx) > 0
                    for j = 1:numel(found_idx)
                        idx = found_idx(j);
                        if db.Age{idx} == age
                            fprintf("Removing %s, %d\n", name, age)
                            db(idx,:) = [];
                            break
                        end
                    end
                end
            end
    end
    save(db_file, "db")
end