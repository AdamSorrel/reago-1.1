-- Created by IntelliJ IDEA.
-- User: stovicek_lab
-- Date: 23/06/17
-- Time: 14:58
-- To change this template use File | Settings | File Templates.
-- KEYS[1] is a list of unique sequences (hashNames)
-- KEYS[2] is a original storage of all equences and keys
-- KEYS[3] is a read position database
-- KEYS[4] is a template position database


redis.replicate_commands()
redis.set_repl(redis.REPL_NONE)

-- cursor is pointing towards the next iteration of HSCAN. When the iterations are finished it returns 0 value again
local cursor = 0

repeat
    -- retrieving cursor from the SCAN method for next iteration
    local output = redis.call("SSCAN", KEYS[1], cursor)

    cursor = output[1]
    if table.getn(output) > 0 then
        local counter = 1
        while counter <= (table.getn(output[2])) do
            local seq = output[2][counter]
            local nameList = redis.call("SMEMBERS", seq)
            local header = table.concat(nameList, '|')

            local readPosition = redis.call("HGET", KEYS[3], nameList[1])
            local templatePosition = redis.call("HGET", KEYS[4], nameList[1])

            redis.call("HSET", KEYS[2], header, seq)
            redis.call("HSET", KEYS[3], header, readPosition)
            redis.call("HSET", KEYS[4], header, templatePosition)

            for position = 1, #nameList do
                -- local otherTemplatePosition = redis.call("HGET", KEYS[4], nameList[position])
                -- if otherTemplatePosition ~= templatePosition then
                --    print(nameList[position])
                --    do return "ERROR: Cannot concatenate - different position occurs in read: " .. nameList[position] end
                -- end
                redis.call("HDEL", KEYS[2], nameList[position])
                redis.call("HDEL", KEYS[3], nameList[position])
                redis.call("HDEL", KEYS[4], nameList[position])
            end
            counter = counter + 1
        end
    end

until tonumber(cursor) == 0

return true

