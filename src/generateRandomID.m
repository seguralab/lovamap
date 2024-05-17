% From ChatGPT to generate a unique identifier

function unique_id = generateRandomID()
    % Generate a unique random identifier
    % Include a high-resolution timestamp to enhance uniqueness
    randNum = rand();  % Generate a random number
    timestamp = char(java.time.Instant.now.toString);  % Get a high-resolution timestamp
    uniqueStr = sprintf('%s-%f', timestamp, randNum);  % Concatenate timestamp and random number
    
    % Use Java's UUID generation for a robust unique ID
    unique_id = char(java.util.UUID.nameUUIDFromBytes(uint8(uniqueStr)));
end
