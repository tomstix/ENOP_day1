% Resetting the database to start fresh
people('peopleDB.mat', 'reset');

% Inserting multiple people into the database
people('peopleDB.mat', 'insert', 'John', 25, 'Anna', 17, 'Mike', 30);

% Listing all entries to verify insertion
people('peopleDB.mat', 'list');

% Removing a specific entry
people('peopleDB.mat', 'remove', 'Anna', 17);

% Listing all entries to verify removal
people('peopleDB.mat', 'list');