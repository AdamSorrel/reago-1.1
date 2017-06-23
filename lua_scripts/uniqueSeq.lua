-- Author: Adam Stovicek
-- Institution: Ben Gurion University of the Negev
-- Date: 20/06/17
-- Time: 19:53
-- To change this template use File | Settings | File Templates.

redis.replicate_commands()
redis.set_repl(redis.REPL_NONE)

-- cursor is pointing towards the next iteration of HSCAN. When the iterations are finished it returns 0 value again
local cursor = 0

repeat
    -- Calling HSCAN on the appropriate database, passing the new cursos value (starting with 0) and approximately
    -- 10 returned value (could be more or less)
    local output = redis.call("HSCAN", KEYS[1], cursor, "COUNT", 10)
    -- Retrieving the cursor of the next iteration
    cursor = output[1]
    -- This is a safety mechanism in case the returned value table is empty
    if table.getn(output) > 0 then
        -- Setting a counter iterator (will iterate only over even numbers)
        local counter = 1
        -- Iteration stops two steps before the end (id ~ end+1 and sequence ~ end+2)
        while (counter <= (table.getn(output[2])-2)) do
            -- Retrieving id of a sequence
            local id = output[2][counter]
            -- Retrieving sequence itself in a +1 position from the local counter
            local seq = output[2][1+counter]

            -- TODO : This might not be necessary
            --if seq == nil then
            --    do return "Sequence nil" end
            --end

            -- Retrieving beginning of the sequence in length specified by ARGV[1]
            -- local subSeq = string.sub(seq, 1, tonumber(KEYS[3]))
            -- TODO : Resolve why is this happening:
            -- Checking if the subset sequence isn't longer then it should be
            --if #subSeq > tonumber(ARGV[1]) then
            --    do return subSeq end
            --end

            -- Saving sequence in a hash with the name of the previously retrieved beginning
            redis.call("SADD", seq, id)
            -- Adding the name of the hash in a set who's name is specified in ARGV[2]
            -- Due to a set property already existing names are not saved in duplicate
            redis.call("SADD", KEYS[2], seq)

            -- Incrementing counter by two position (skipping sequence position)
            counter = counter + 2
        end
    end

until tonumber(cursor) == 0

return true