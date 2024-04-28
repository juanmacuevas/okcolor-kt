package nl.juanma

import java.lang.Math.pow

class OkColor {


    fun okHslToSrgb(h: Double, s: Double, l: Double): Triple<Int, Int, Int> {
        if (l == 1.0) {
            return Triple(255, 255, 255)
        } else if (l == 0.0) {
            return Triple(0, 0, 0)
        }

        val a = Math.cos(2 * Math.PI * h)
        val b = Math.sin(2 * Math.PI * h)
        val L = toeInv(l)

        val Cs = getCs(L, a, b)
        val C0 = Cs[0]
        val Cmid = Cs[1]
        val Cmax = Cs[2]

        val (C, t, k0, k1, k2) = if (s < 0.8) {
            val t = 1.25 * s
            val k0 = 0.0
            val k1 = 0.8 * C0
            val k2 = 1 - k1 / Cmid
            listOf(k0 + t * k1 / (1 - k2 * t), t, k0, k1, k2)
        } else {
            val t = 5 * (s - 0.8)
            val k0 = Cmid
            val k1 = 0.2 * Cmid * Cmid * 1.25 * 1.25 / C0
            val k2 = 1 - (k1) / (Cmax - Cmid)
            listOf(k0 + t * k1 / (1 - k2 * t), t, k0, k1, k2)
        }

        val rgb = oklabToLinearSrgb(L, C * a, C * b)
        return Triple(
            (255 * srgbTransferFunction(rgb.first)).toInt(),
            (255 * srgbTransferFunction(rgb.second)).toInt(),
            (255 * srgbTransferFunction(rgb.third)).toInt()
        )
    }

    fun oklabToLinearSrgb(L: Double, a: Double, b: Double): Triple<Double, Double, Double> {
        val l_ = L + 0.3963377774 * a + 0.2158037573 * b
        val m_ = L - 0.1055613458 * a - 0.0638541728 * b
        val s_ = L - 0.0894841775 * a - 1.2914855480 * b

        val l = l_ * l_ * l_
        val m = m_ * m_ * m_
        val s = s_ * s_ * s_

        return Triple(
            4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s,
            -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
            -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s
        )
    }

    internal fun toe(x: Double): Double {
        val k1 = 0.206
        val k2 = 0.03
        val k3 = (1 + k1) / (1 + k2)

        return 0.5 * (k3 * x - k1 + Math.sqrt(
            (k3 * x - k1) * (k3 * x - k1) + 4 * k2 * k3 * x
        ))
    }

    internal fun toeInv(x: Double): Double {
        val k1 = 0.206
        val k2 = 0.03
        val k3 = (1 + k1) / (1 + k2)

        return (x * x + k1 * x) / (k3 * (x + k2))
    }

    // Finds the maximum saturation possible for a given hue that fits in sRGB
    // Saturation here is defined as S = C/L
    //  a and b must be normalized so a^2 + b^2 == 1
    internal fun computeMaxSaturation(a: Double, b: Double): Double {
        // Max saturation will be when one of r, g, or b goes below zero.
        // Select different coefficients depending on which component goes below zero first
        var k0: Double
        var k1: Double
        var k2: Double
        var k3: Double
        var k4: Double
        var wl: Double
        var wm: Double
        var ws: Double

        if (-1.88170328 * a - 0.80936493 * b > 1) {
            // Red component
            k0 = 1.19086277; k1 = 1.76576728; k2 = 0.59662641; k3 = 0.75515197; k4 = 0.56771245
            wl = 4.0767416621; wm = -3.3077115913; ws = 0.2309699292
        } else if (1.81444104 * a - 1.19445276 * b > 1) {
            // Green component
            k0 = 0.73956515; k1 = -0.45954404; k2 = 0.08285427; k3 = 0.12541070; k4 = 0.14503204
            wl = -1.2684380046; wm = 2.6097574011; ws = -0.3413193965
        } else {
            // Blue component
            k0 = 1.35733652; k1 = -0.00915799; k2 = -1.15130210; k3 = -0.50559606; k4 = 0.00692167
            wl = -0.0041960863; wm = -0.7034186147; ws = 1.7076147010
        }

        // Approximate max saturation using a polynomial:
        var S = k0 + k1 * a + k2 * b + k3 * a * a + k4 * a * b

        // Do one step Halley's method to get closer
        val kL = 0.3963377774 * a + 0.2158037573 * b
        val kM = -0.1055613458 * a - 0.0638541728 * b
        val kS = -0.0894841775 * a - 1.2914855480 * b

        val l_ = 1 + S * kL
        val m_ = 1 + S * kM
        val s_ = 1 + S * kS

        val l = l_ * l_ * l_
        val m = m_ * m_ * m_
        val s = s_ * s_ * s_

        val lDS = 3 * kL * l_ * l_
        val mDS = 3 * kM * m_ * m_
        val sDS = 3 * kS * s_ * s_

        val lDS2 = 6 * kL * kL * l_
        val mDS2 = 6 * kM * kM * m_
        val sDS2 = 6 * kS * kS * s_

        val f = wl * l + wm * m + ws * s
        val f1 = wl * lDS + wm * mDS + ws * sDS
        val f2 = wl * lDS2 + wm * mDS2 + ws * sDS2

        S -= f * f1 / (f1 * f1 - 0.5 * f * f2)

        return S
    }



    internal fun findCusp(a: Double, b: Double): Pair<Double, Double> {
        // First, find the maximum saturation (saturation S = C/L)
        val S_cusp = computeMaxSaturation(a, b)

        // Convert to linear sRGB to find the first point where at least one of r, g, or b >= 1:
        val rgbAtMax = oklabToLinearSrgb(1.0, S_cusp * a, S_cusp * b)
        val L_cusp = Math.cbrt(1.0 / maxOf(rgbAtMax.first, rgbAtMax.second, rgbAtMax.third))
        val C_cusp = L_cusp * S_cusp

        return Pair(L_cusp, C_cusp)
    }

    internal fun findGamutIntersection(
        a: Double, b: Double, L1: Double, C1: Double, L0: Double, cusp: Pair<Double, Double>? = null
    ): Double {
        val cuspEffective = cusp ?: findCusp(a, b)

        val (L_cusp, C_cusp) = cuspEffective
        var t: Double

        if ((L1 - L0) * C_cusp - (L_cusp - L0) * C1 <= 0) {
            // Lower half
            t = C_cusp * L0 / (C1 * L_cusp + C_cusp * (L0 - L1))
        } else {
            // Upper half
            // First intersect with triangle
            t = C_cusp * (L0 - 1) / (C1 * (L_cusp - 1) + C_cusp * (L0 - L1))

            // Halley's method iteration
            val dL = L1 - L0
            val dC = C1

            val kL = 0.3963377774 * a + 0.2158037573 * b
            val kM = -0.1055613458 * a - 0.0638541728 * b
            val kS = -0.0894841775 * a - 1.2914855480 * b

            val lDt = dL + dC * kL
            val mDt = dL + dC * kM
            val sDt = dL + dC * kS

            var L = L0 * (1 - t) + t * L1
            var C = t * C1

            val l_ = L + C * kL
            val m_ = L + C * kM
            val s_ = L + C * kS

            val l = l_ * l_ * l_
            val m = m_ * m_ * m_
            val s = s_ * s_ * s_

            val ldt = 3 * lDt * l_ * l_
            val mdt = 3 * mDt * m_ * m_
            val sdt = 3 * sDt * s_ * s_

            val ldt2 = 6 * lDt * lDt * l_
            val mdt2 = 6 * mDt * mDt * m_
            val sdt2 = 6 * sDt * sDt * s_

            var r = 4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s - 1
            var r1 = 4.0767416621 * ldt - 3.3077115913 * mdt + 0.2309699292 * sdt
            var r2 = 4.0767416621 * ldt2 - 3.3077115913 * mdt2 + 0.2309699292 * sdt2

            var u_r = r1 / (r1 * r1 - 0.5 * r * r2)
            var t_r = -r * u_r

            t += kotlin.math.min(t_r, 10e5) // Adjust this line based on which color component (r, g, b) is relevant
        }

        return t
    }

    internal fun getSTMax(a: Double, b: Double, cusp: Pair<Double, Double>? = null): Pair<Double, Double> {
        val effectiveCusp = cusp ?: findCusp(a, b)
        val (L, C) = effectiveCusp
        return Pair(C / L, C / (1 - L))
    }


    internal fun getSTMid(a: Double, b: Double): Pair<Double, Double> {
        val S =
            0.11516993 + 1 / (7.44778970 + 4.15901240 * b + a * (-2.19557347 + 1.75198401 * b + a * (-2.13704948 - 10.02301043 * b + a * (-4.24894561 + 5.38770819 * b + 4.69891013 * a))))

        val T =
            0.11239642 + 1 / (1.61320320 - 0.68124379 * b + a * (0.40370612 + 0.90148123 * b + a * (-0.27087943 + 0.61223990 * b + a * (0.00299215 - 0.45399568 * b - 0.14661872 * a))))

        return Pair(S, T)
    }

    internal fun srgbTransferFunction(a: Double): Double {
        return if (0.0031308 >= a) 12.92 * a else 1.055 * pow(a, .4166666666666667) - .055
    }

    internal fun srgbTransferFunction_inv(a: Double): Double {
        return if (.04045 < a) pow((a + .055) / 1.055, 2.4) else a / 12.92
    }



    internal fun getCs(L: Double, a: Double, b: Double): List<Double> {
        val cusp = findCusp(a, b)
        val Cmax = findGamutIntersection(a, b, L, 1.0, L, cusp)
        val STmax = getSTMax(a, b, cusp)

        val Smid =
            0.11516993 + 1 / (7.44778970 + 4.15901240 * b + a * (-2.19557347 + 1.75198401 * b + a * (-2.13704948 - 10.02301043 * b + a * (-4.24894561 + 5.38770819 * b + 4.69891013 * a))))

        val Tmid =
            0.11239642 + 1 / (1.61320320 - 0.68124379 * b + a * (0.40370612 + 0.90148123 * b + a * (-0.27087943 + 0.61223990 * b + a * (0.00299215 - 0.45399568 * b - 0.14661872 * a))))

        val k = Cmax / minOf((L * STmax.first), (1 - L) * STmax.second)

        val Ca = L * Smid
        val Cb = (1 - L) * Tmid
        val Cmid = 0.9 * k * Math.sqrt(Math.sqrt(1 / (1 / (Ca * Ca * Ca * Ca) + 1 / (Cb * Cb * Cb * Cb))))

        val Ca0 = L * 0.4
        val Cb0 = (1 - L) * 0.8
        val C0 = Math.sqrt(1 / (1 / (Ca0 * Ca0) + 1 / (Cb0 * Cb0)))

        return listOf(C0, Cmid, Cmax)
    }


}