-- Created by IntelliJ IDEA.
-- User: stovicek_lab
-- Date: 23/06/17
-- Time: 14:58
-- To change this template use File | Settings | File Templates.
-- KEYS[1] is a list of unique sequences (hashNames)
-- KEYS[2] is a original storage of all equences and keys (read_sequence)

-- KEYS[3] is a read position database (read_position)
-- KEYS[4] is a template position database (template_position)

-- Local function for splitting strings using a 'split' character
local function split(inputstr, sep)
        if sep == nil then
                sep = "%s"
        end
        local t={} ; local i=1
        for str in string.gmatch(inputstr, "([^"..sep.."]+)") do
                t[i] = str
                i = i + 1
        end
        return t
end

redis.replicate_commands()
redis.set_repl(redis.REPL_NONE)

-- cursor is pointing towards the next iteration of HSCAN. When the iterations are finished it returns 0 value again
local cursor = 0

repeat
    -- retrieving cursor from the SCAN method for next iteration
    local output = redis.call("SSCAN", KEYS[1], cursor)
    -- saving a cursor for the next iteration
    cursor = output[1]
    -- This is just a safety mechanism, should the returned list be empty
    if table.getn(output) > 0 then

        -- counter to iterate across the returned table
        local counter = 1
        while counter <= (table.getn(output[2])) do
            -- The object output is nested with cursor at position 1 and a table with sequence list at position 2
            local seq = output[2][counter]
            -- retrieving a list of headers from REDIS set
            local headerList = redis.call("SMEMBERS", seq)

            -- defining persistent variables for the next loop
            local nameList = {}
            local headerParameters = ''

            -- looping over the headerList
            for seqNum = 1, #headerList do
                -- splitting sequence header at ';'
                local splitSeq = split(headerList[seqNum], ';')
                -- appending a new name to the nameList
                table.insert(nameList, splitSeq[1])
                -- retrieving the rest of parameters (assumed identical for all headers)
                headerParameters = ';' .. splitSeq[2] .. ';' .. splitSeq[3] .. ';' .. splitSeq[4]
            end

            -- concatenating header names separated with the '|'
            local headerNames = table.concat(nameList, '|')
            -- concatenating the header names with the last headerParaemters (assumed all identical)
            local header = headerNames .. headerParameters

            --local readPosition = redis.call("HGET", KEYS[3], nameList[1])
            --local templatePosition = redis.call("HGET", KEYS[4], nameList[1])

            redis.call("HSET", KEYS[2], header, seq)
            --redis.call("HSET", KEYS[3], header, readPosition)
            --redis.call("HSET", KEYS[4], header, templatePosition)

            -- cleaning the database
            for position = 1, #nameList do
                -- local otherTemplatePosition = redis.call("HGET", KEYS[4], nameList[position])
                -- if otherTemplatePosition ~= templatePosition then
                --    print(nameList[position])
                --    do return "ERROR: Cannot concatenate - different position occurs in read: " .. nameList[position] end
                -- end

                redis.call("HDEL", KEYS[2], nameList[position])


                --redis.call("HDEL", KEYS[3], nameList[position])
                --redis.call("HDEL", KEYS[4], nameList[position])
            end
            counter = counter + 1
        end
    end

until tonumber(cursor) == 0

return true

