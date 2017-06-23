--
-- Created by IntelliJ IDEA.
-- User: stovicek_lab
-- Date: 23/06/17
-- Time: 14:58
-- To change this template use File | Settings | File Templates.
--

redis.replicate_commands()
redis.set_repl(redis.REPL_NONE)

-- cursor is pointing towards the next iteration of HSCAN. When the iterations are finished it returns 0 value again
local cursor = 0

repeat


until tonumber(cursor) == 0

return true

