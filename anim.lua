-------------------------------------------------------------------------
--
-- A module for tweening Roblox datatypes.
--
--
-- Copyright 2017 Parker Stebbins <parker@fractality.io>
--
-- Permission is hereby granted, free of charge, to any person obtaining a copy of this software
-- and associated documentation files (the "Software"), to deal in the Software without
-- restriction, including without limitation the rights to use, copy, modify, merge, publish,
-- distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the
-- Software is furnished to do so, subject to the following conditions:
--
-- The above copyright notice and this permission notice shall be included in all copies or
-- substantial portions of the Software.
--
-- THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
-- BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
-- NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
-- DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
-- OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
--
--
-- API ------------------------------------------------------------------
--
-- Animate(Instance obj, dictionary<string, Variant> plist, double time, string style, int priority = 0)
-- | @brief Animates properties of obj to the given values.
-- | @param obj       The object to animate.
-- | @param plist     Dictionary of target properties.
-- | @param time      The duration of the animation.
-- | @param style     Name of an easing style.
-- | @param priority  The priority of the animation. New animations will only interrupt running animations of lower or equal priority.
--
-- AnimateAsync(Instance obj, dictionary<string, Variant> plist, double time, string style, int priority=0[, function callback()])
-- | @brief Non-blocking equivalent of Animate.
-- | @param callback  A callback function to run after the animation stops.
--
-- bool Interrupt(Instance obj[, int priority])
-- | @brief Stop the animation (if any) acting on this object.
-- | @param obj  The object whose animation will be interrupted.
-- | @param priority  If supplied, only animations under or equal to priority will be interrupted.
--
--
-- Easing styles --------------------------------------------------------
--
-- | linear    |
-- | inQuad    | outQuad    | inOutQuad    |
-- | inCubic   | outCubic   | inOutCubic   |
-- | inQuart   | outQuart   | inOutQuart   |
-- | inQuint   | outQuint   | inOutQuint   |
-- | inElastic | outElastic | inOutElastic |
-- | inBack    | outBack    | inOutBack    |
-- | inBounce  | outBounce  |     n/a      |
--
-- in:    start slow and speed up
-- out:   start fast and slow down
-- inOut: start slow, speed up, then slow down
--
--[[ Usage example ------------------------------------------------------

anim.Animate(frame, {
	Position = UDim2.new(0.5, 0, 0.5, 0);
	Size     = UDim2.new(0, 0, 0, 0);
})

anim.AnimateAsync(frame, {
	Position = UDim2.new(0.3, 0, 0.3, 0);
	Size     = UDim2.new(0.4, 0, 0.4, 0);
}, 0.5, 'outBack')

wait(0.25)

anim.Animate(frame, {
	Position = UDim2.new(0.3, 0, 0.3, 0);
	Size     = UDim2.new(0.4, 0, 0.4, 0);
}, 0.5, 'outBack')

-----------------------------------------------------------------------]]

local DEFAULT_PRIORITY = 0
local TYPE_ASSERTS = true

local api = {}
local easingStyles = {}

do
	local typeDispatch = setmetatable({}, {
		__index = function(t, k)
			error(('No interpolator is defined for type %s.'):format(tostring(k)), 4)
		end
	})
	do
		-- bool
		function typeDispatch.bool(v0, v1)
			return function(t)
				if t < 0.5 then
					return v0
				else
					return v1
				end
			end
		end

		-- number
		function typeDispatch.number(v0, v1)
			local dv = v1 - v0
			return function(t)
				return v0 + dv*t
			end
		end

		do -- string
			local s_match = string.match
			local s_format = string.format
			local atof = tonumber

			function typeDispatch.string(v0, v1)
				local n0, d do
					local sign0, h0, m0, s0 = s_match(v0, '^([+-]?)(%d*):[+-]?(%d*):[+-]?(%d*)$')
					local sign1, h1, m1, s1 = s_match(v1, '^([+-]?)(%d*):[+-]?(%d*):[+-]?(%d*)$')
					if sign0 and sign1 then
						n0       = 3600*(atof(h0) or 0) + 60*(atof(m0) or 0) + (atof(s0) or 0)
						local n1 = 3600*(atof(h1) or 0) + 60*(atof(m1) or 0) + (atof(s1) or 0)
						if sign0 == '-' then
							n0 = -n0
						end
						d = (43200 + (sign1 ~= '-' and n1 or -n1) - n0)%86400 - 43200
					else
						error('Invalid TimeOfDay string', 4)
					end
				end

				return function(t)
					local fs = (n0 + d*t)%86400
					local s = fs > 0 and fs or -fs
					return s_format(
						fs < 0 and '-%.2u:%.2u:%.2u' or '%.2u:%.2u:%.2u',
						(s - s%3600)/3600,
						(s%3600 - s%60)/60,
						s%60
					)
				end
			end
		end

		do -- CFrame
			local Slerp = CFrame.new().lerp
			function typeDispatch.CFrame(v0, v1)
				return function(t)
					return Slerp(v0, v1, t)
				end
			end
		end

		do -- Color3
			local clamp = math.clamp
			local C3 = Color3.new
			local black = C3(0, 0, 0)

			local function RgbToLuv13(c)
				local r, g, b = c.r, c.g, c.b
				r = r < 0.0404482362771076 and r/12.92 or 0.87941546140213*(r + 0.055)^2.4
				g = g < 0.0404482362771076 and g/12.92 or 0.87941546140213*(g + 0.055)^2.4
				b = b < 0.0404482362771076 and b/12.92 or 0.87941546140213*(b + 0.055)^2.4
				local y = 0.2125862307855956*r + 0.71517030370341085*g + 0.0722004986433362*b
				local z = 3.6590806972265883*r + 11.4426895800574232*g + 4.1149915024264843*b
				local l = y > 0.008856451679035631 and 116*y^(1/3) - 16 or 903.296296296296*y
				if z > 1e-15 then
					local x = 0.9257063972951867*r - 0.8333736323779866*g - 0.09209820666085898*b
					return l, l*x/z, l*(9*y/z - 0.46832)
				else
					return l, -0.19783*l, -0.46832*l
				end
			end

			function typeDispatch.Color3(c0, c1)
				local l0, u0, v0 = RgbToLuv13(c0)
				local l1, u1, v1 = RgbToLuv13(c1)

				c0, c1 = nil, nil

				return function(t)
					local l = (1 - t)*l0 + t*l1
					if l < 0.0197955 then
						return black
					end

					local u = ((1 - t)*u0 + t*u1)/l + 0.19783
					local v = ((1 - t)*v0 + t*v1)/l + 0.46832

					local y = (l + 16)/116
					y = y > 0.206896551724137931 and y*y*y or 0.12841854934601665*y - 0.01771290335807126
					local x = y*u/v
					local z = y*((3 - 0.75*u)/v - 5)

					local r =  7.2914074*x - 1.5372080*y - 0.4986286*z
					local g = -2.1800940*x + 1.8757561*y + 0.0415175*z
					local b =  0.1253477*x - 0.2040211*y + 1.0569959*z

					if r < 0 and r < g and r < b then
						r, g, b = 0, g - r, b - r
					elseif g < 0 and g < b then
						r, g, b = r - g, 0, b - g
					elseif b < 0 then
						r, g, b = r - b, g - b, 0
					end

					return C3(
						clamp(r < 3.1306684425e-3 and 12.92*r or 1.055*r^(1/2.4) - 0.055, 0, 1),
						clamp(g < 3.1306684425e-3 and 12.92*g or 1.055*g^(1/2.4) - 0.055, 0, 1),
						clamp(b < 3.1306684425e-3 and 12.92*b or 1.055*b^(1/2.4) - 0.055, 0, 1)
					)
				end
			end
		end


		do -- NumberRange
			local NR = NumberRange.new
			function typeDispatch.NumberRange(v0, v1)
				local min0, max0 = v0.Min, v0.Max
				local dmin, dmax = v1.Min - min0, v1.Max - max0
				v0, v1 = nil, nil
				return function(t)
					return NR(min0 + t*dmin, max0 + t*dmax)
				end
			end
		end

		do -- NumberSequenceKeypoint
			local NSK = NumberSequenceKeypoint.new
			function typeDispatch.NumberSequenceKeypoint(v0, v1)
				local t0, v0, e0 = v0.Time, v0.Value, v0.Envelope
				local dt, dv, de = v1.Time - t0, v1.Value - v0, v1.Envelope - e0
				v1 = nil
				return function(t)
					return NSK(t0 + t*dt, v0 + t*dv, e0 + t*de)
				end
			end
		end

		do -- PhysicalProperties
			local PP = PhysicalProperties.new
			function typeDispatch.PhysicalProperties(v0, v1)
				local d0, e0, ew0, f0, fw0 =
					v0.Density,
					v0.Elasticity,
					v0.ElasticityWeight,
					v0.Friction,
					v0.FrictionWeight
				local dd, de, dew, df, dfw =
					v1.Density - d0,
					v1.Elasticity - e0,
					v1.ElasticityWeight - ew0,
					v1.Friction - f0,
					v1.FrictionWeight - fw0
				v0, v1 = nil, nil
				return function(t)
					return PP(d0 + t*dd, e0 + t*de, ew0 + t*dew, f0 + t*df, fw0 + t*dfw)
				end
			end
		end

		do -- Ray
			local R = Ray.new
			local V3 = Vector3.new
			function typeDispatch.Ray(v0, v1)
				local o0, d0, o1, d1 =
					v0.Origin, v0.Direction,
					v1.Origin, v1.Direction
				local ox0, oy0, oz0, dx0, dy0, dz0, dx1, dy1, dz1 =
					o0.x, o0.y, o0.z,
					d0.x, d0.y, d0.z,
					d1.x, d1.y, d1.z
				local dox, doy, doz, ddx, ddy, ddz =
					o1.x - ox0, o1.y - oy0, o1.z - oz0,
					d1.x - dx0, d1.y - dy0, d1.z - dz0
				v0, v1, o0, d0, o1, d1 = nil, nil, nil, nil, nil, nil
				return function(t)
					return R(
						V3(ox0 + t*dox, oy0 + t*doy, oz0 + t*doz),
						V3(dx0 + t*ddx, dy0 + t*ddy, dz0 + t*ddz)
					)
				end
			end
		end

		do -- UDim
			local UD = UDim.new
			function typeDispatch.UDim(v0, v1)
				local sc, of = v0.Scale, v0.Offset
				local dsc, dof = v1.Scale - sc, v1.Offset - of
				v0, v1 = nil, nil
				return function(t)
					return UD(sc + t*dsc, of + t*dof)
				end
			end
		end

		do -- UDim2
			local Lerp = UDim2.new().Lerp
			function typeDispatch.UDim2(v0, v1)
				return function(t)
					return Lerp(v0, v1, t)
				end
			end
		end

		do -- Vector2
			local V2 = Vector2.new
			function typeDispatch.Vector2(v0, v1)
				local x, y = v0.x, v0.y
				local dx, dy = v1.x - x, v1.y - y
				v0, v1 = nil, nil
				return function(t)
					return V2(x + t*dx, y + t*dy)
				end
			end
		end

		do -- Vector3
			local V3 = Vector3.new
			function typeDispatch.Vector3(v0, v1)
				local x, y, z = v0.x, v0.y, v0.z
				local dx, dy, dz = v1.x - x, v1.y - y, v1.z - z
				v0, v1 = nil, nil
				return function(t)
					return V3(x + t*dx, y + t*dy, z + t*dz)
				end
			end
		end
	end

	local vlock = {} -- lock record
	local plock = {} -- priority record

	local Interpolate do
		local RunService = game:GetService'RunService'
		local step
		if RunService:IsServer() and not RunService:IsClient() then
			step = RunService.Heartbeat
		else
			step = RunService.RenderStepped
		end

		local Await = step.Wait
		local next = next
		local tick = tick
		local typeof = typeof

		local dtRender = 1/60 do
			local last = tick()
			step:Connect(function()
				local t = tick()
				last, dtRender = t, dtRender + (t - last - dtRender)/4
			end)
		end

		function Interpolate(obj, props, t, Ease, priority, Callback)
			local priMap = plock[obj]

			if not t or t < dtRender then
				for prop, val in next, props do
					local overridenPriority = priMap and priMap[prop]
					if (overridenPriority or -1/0) <= priority then
						local excludeMap = vlock[prop]
						if excludeMap then
							local mv = excludeMap[obj]
							if mv then
								excludeMap[obj] = (mv + 1)%2^53
							end
						end
						if overridenPriority then
							priMap[prop] = nil
						end
						obj[prop] = val
					end
				end
				if Callback then
					return Callback()
				end
				return
			elseif type(Ease) ~= 'function' then
				Ease = easingStyles[Ease]
			end

			local t0 = tick() - dtRender
			local x = Ease(dtRender/t)
			local flerps = {}
			local mxkeys = {}
			local noscan = true

			if not priMap then
				priMap = {}
				noscan = false
				plock[obj] = priMap
			end

			for prop, val in next, props do
				if noscan or (priMap[prop] or -1/0) <= priority then
					local excludeMap = vlock[prop]
					if not excludeMap then
						excludeMap = {}
						vlock[prop] = excludeMap
					end
					local mv = ((excludeMap[obj] or 0) + 1)%2^53
					local typeLerp = typeDispatch[typeof(val)](obj[prop], val)
					mxkeys[prop] = mv
					excludeMap[obj] = mv
					priMap[prop] = priority
					flerps[prop] = typeLerp
					obj[prop] = typeLerp(x)
				end
			end

			repeat
				Await(step)
				local elapsed = tick() - t0
				if elapsed >= t then
					break
				end
				x = Ease(elapsed/t)
				for prop, Flerp in next, flerps do
					if vlock[prop][obj] ~= mxkeys[prop] then
						flerps[prop] = nil
					else
						obj[prop] = Flerp(x)
					end
				end
			until not next(flerps)

			for prop, val in next, flerps do
				local mx = vlock[prop]
				obj[prop] = props[prop]
				mx[obj] = nil
				priMap[prop] = nil
			end
			if not next(priMap) then
				plock[obj] = nil
			end
			if Callback then
				return Callback()
			end
		end
	end

	local type = type
	local c_wrap = coroutine.wrap

	local chktype
	if TYPE_ASSERTS then
		local select = select

		local function GetCallerName(lvl)
			local f = debug.traceback():gmatch'Script \'.-\', Line %d+ %- [globacupvefidmth]+ ([_%a][_%w]*)\n'
			for _ = 1, lvl or 1 do
				f()
			end
			return f()
		end

		function chktype(argno, dat, ...)
			local tl = {...}
			local td = type(dat)
			local er = true
			for i = 1, select('#', ...) do
				er = er and (tl[i] ~= td)
			end
			if er then
				error(('bad argument #%u to %s (%s expected, got %s)'):format(argno, GetCallerName(2), table.concat(tl, ' or '), td), 3)
			end
		end
	end

	function api.Animate(obj, props, time, style, priority)
		if TYPE_ASSERTS then
			chktype(1, obj,      'userdata', 'table')
			chktype(2, props,    'table')
			chktype(3, time,     'nil', 'number')
			chktype(4, style,    'nil', 'string', 'function')
			chktype(5, priority, 'nil', 'number')
		end
		return Interpolate(obj, props, time, style, priority or DEFAULT_PRIORITY)
	end

	function api.AnimateAsync(obj, props, time, style, priority, Callback)
		if TYPE_ASSERTS then
			chktype(1, obj,      'userdata', 'table')
			chktype(2, props,    'table')
			chktype(3, time,     'nil', 'number')
			chktype(4, style,    'nil', 'string', 'function')
			chktype(5, priority, 'nil', 'number')
			chktype(6, Callback, 'nil', 'function')
		end
		return c_wrap(Interpolate)(obj, props, time, style, priority or DEFAULT_PRIORITY, Callback)
	end

	function api.Interrupt(obj)
		if TYPE_ASSERTS then
			chktype(1, obj, 'nil', 'userdata', 'table')
		end
		for _, s in next, vlock do
			if not obj or s[obj] then
				s[obj] = nil
			end
		end
		if obj then
			plock[obj] = nil
		else
			plock = {}
		end
	end
end


do
	local cos = math.cos
	local sin = math.sin
	local pi = math.pi
	local psi = pi/2
	local tau = pi*2
	local exp = math.exp

	function easingStyles.linear(t)
		return t
	end

	do -- Quad
		function easingStyles.inQuad(t)
			return t*t
		end

		function easingStyles.outQuad(t)
			return t*(2 - t)
		end

		function easingStyles.inOutQuad(t)
			if t < 0.5 then
				return 2*t*t
			else
				return 2*(2 - t)*t - 1
			end
		end
	end

	do -- Cubic
		function easingStyles.inCube(t)
			return t*t*t
		end

		function easingStyles.outCube(t)
			return 1 - (1 - t)*(1 - t)*(1 - t)
		end

		function easingStyles.inOutCube(t)
			if t < 0.5 then
				return 4*t*t*t
			else
				return 1 + 4*(t - 1)*(t - 1)*(t - 1)
			end
		end
	end

	do -- Quart
		function easingStyles.inQuart(t)
			return (t*t)*(t*t)
		end

		function easingStyles.outQuart(t)
			return 1 - ((t - 1)*(t - 1))*((t - 1)*(t - 1))
		end

		function easingStyles.inOutQuart(t)
			if t < 0.5 then
				return 8*(t*t)*(t*t)
			else
				return 1 - 8*((t - 1)*(t - 1))*((t - 1)*(t - 1))
			end
		end
	end

	do -- Quint
		function easingStyles.inQuint(t)
			return (t*t)*(t*t)*t
		end

		function easingStyles.outQuint(t)
			t = t - 1
			return (t*t)*(t*t)*t + 1
		end

		function easingStyles.inOutQuint(t)
			if t < 0.5 then
				return 16*(t*t)*(t*t)*t
			else
				t = t - 1
				return 16*(t*t)*(t*t)*t + 1
			end
		end
	end

	do -- Back
		local o = 2 -- Overshoot

		function easingStyles.inBack(t)
			return t*t*((o + 1)*t - o)
		end

		function easingStyles.outBack(t)
			return (t - 1)*(t - 1)*(t*o + t - 1) + 1
		end

		function easingStyles.inOutBack(t)
			if t < 0.5 then
				return 2*t*t*(2*(1 + o)*t - o)
			else
				return 1 + 2*(t - 1)*(t - 1)*(2*(1 + o)*t - 2 - o)
			end
		end
	end

	do -- Sine
		function easingStyles.inSine(t)
			return 1 - cos(t*psi)
		end

		function easingStyles.outSine(t)
			return sin(t*psi)
		end

		function easingStyles.inOutSine(t)
			return (1 - cos(pi*t))/2
		end
	end

	do -- Bounce
		function easingStyles.outBounce(t)
			if t < 4/11 then
				return (121/16)*t*t
			elseif t < 8/11 then
				return 3 + t*(11*t - 12)*(11/16)
			elseif t < 10/11 then
				return 6 + t*(11*t - 18)*(11/16)
			else
				return 63/8 + t*(11*t - 21)*(11/16)
			end
		end

		function easingStyles.inBounce(t)
			if t > 7/11 then
				return 1 - (t - 1)*(t - 1)*(121/16)
			elseif t > 3/11 then
				return (11*t - 7)*(11*t - 3)/-16
			elseif t > 1/11 then
				return (11*(4 - 11*t)*t - 3)/16
			else
				return t*(11*t - 1)*(-11/16)
			end
		end
	end

	do -- Elastic
		local g = 8 -- Amplitude falloff
		local h = 4 -- Number of peaks

		-- Calculate a time multiplier so that the derivative at the end of the wobble is zero.
		-- Then, scale the output so that f(0)=1.
		local hScale, vScale do
			h = pi*(2*h - 3/2)
			local i = 0
			local prev = 0
			local x = 1
			repeat
				local s, c = sin(h*x), cos(h*x)
				prev, x = x, x + (2*h*(1 + x*g)*c + (g*(2 + x*g) - h*h*x)*s)/(h*(h*h*x - 3*g*(2 + x*g))*c + (3*h*h*(1 + x*g) - g*g*(3 + x*g))*s)
				i = i + 1
			until x >= prev or i > 63
			hScale, vScale = x, 1/(exp((x - 1)*g)*x*sin(h*x))
		end

		function easingStyles.inElastic(t)
			return exp((t*hScale - 1)*g)*t*hScale*sin(h*t*hScale)*vScale
		end

		function easingStyles.outElastic(t)
			return 1 + (exp(g*(hScale - hScale*t - 1))*hScale*(t - 1)*sin(h*hScale*(1 - t)))*vScale
		end

		function easingStyles.inOutElastic(t)
			if t < 0.5 then
				return (exp(g*(2*hScale*t - 1))*hScale*t*sin(2*h*hScale*t))*vScale
			else
				return 1 + (exp(g*(hScale*(2 - 2*t) - 1))*hScale*(t - 1)*sin(h*hScale*(2 - 2*t)))*vScale
			end
		end
	end

	do -- Expo
		local s = 4

		function easingStyles.inExpo(t)
			return t*t*exp(s*(t - 1))
		end

		function easingStyles.outExpo(t)
			return 1 - (1 - t)*(1 - t)/exp(s*t)
		end

		function easingStyles.inOutExpo(t)
			if t < 0.5 then
				return 2*t*t*exp(s*(2*t - 1))
			else
				return 1 - 2*(t - 1)*(t - 1)*exp(s*(1 - 2*t))
			end
		end
	end

	setmetatable(easingStyles, {
		__index = function(t, k)
			error(tostring(k) .. ' is not a valid easing style.', 2)
		end;
	})
end

return api
