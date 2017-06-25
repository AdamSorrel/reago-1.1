--
-- Created by IntelliJ IDEA.
-- User: stovicek_lab
-- Date: 23/06/17
-- Time: 14:58
-- To change this template use File | Settings | File Templates.
--

-- KEYS[1] is a list of unique sequences
-- KEYS[2] is a original storage of all equences and keys


redis.replicate_commands()
redis.set_repl(redis.REPL_NONE)

-- cursor is pointing towards the next iteration of HSCAN. When the iterations are finished it returns 0 value again
local cursor = 0

repeat
    local output = redis.call("SSCAN", KEYS[1], cursor)
    cursor = output[1]
    if table.getn(output) > 0 then
        local counter = 1
        while counter <= (table.getn(output[2])) do
            local seq = output[2][counter]
            local nameList = redis.call("SMEMBERS", seq)
            local header = table.concat(nameList, '|')

            redis.call("HSET", KEYS[2] ,header, seq)
            for position = 1, #nameList do
                redis.call("HDEL", KEYS[2], nameList[position])
            end
            counter = counter + 1
        end
    end

until tonumber(cursor) == 0

return true

