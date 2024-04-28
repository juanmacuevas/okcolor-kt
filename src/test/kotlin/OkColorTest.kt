package nl.juanma

import nl.juanma.OkColor
import org.junit.jupiter.api.Assertions.*
import kotlin.math.abs
import kotlin.test.Test

class OkColorTest {

    @Test
    fun okhslToRgb() {

    }

    @Test
    fun rgbToOkhsl() {
    }



    @Test
    fun testToe() {
        val x = 10.0
        val expected = 11.533 // Replace with your expected value
        val result = OkColor().toe(x)
        assertEquals(expected, result, 0.01) // Allow for a small error

    }

    @Test
    fun testToeInv() {
        val x = 5.0
        val expected = 4.419 // Replace with your expected value
        val result = OkColor().toeInv(x)
        assertEquals(expected, result, 0.01) // Allow for a small error margin
    }

    @Test
    fun testComputeMaxSaturation() {
        // Test data
        val testCases = arrayOf(
            Pair(1.0 / Math.sqrt(2.0), 1.0 / Math.sqrt(2.0)), // Example where a^2 + b^2 == 1
            Pair(0.5, 0.5),                                    // Non-normalized values
            Pair(0.0, 1.0),                                    // Edge case on the unit circle
            Pair(1.0, 0.0)                                     // Edge case on the unit circle
        )

        // Expected results (These should be filled in with expected values for testing)
        val expectedResults = arrayOf(
            0.28590006172378446,
            0.42172089988292466,
            0.204357271698336,
            0.4053912761182488
        )

        // Test loop
        for ((index, testCase) in testCases.withIndex()) {
            val (a, b) = testCase
            val result = OkColor().computeMaxSaturation(a, b)
            assert(result == expectedResults[index]) {
                "Test failed for input (a: $a, b: $b). Expected ${expectedResults[index]}, but got $result"
            }
            println("Test passed for input (a: $a, b: $b). Expected and actual value: $result")
        }
    }

    @Test
    fun testOklabToLinearSrgb() {
        val L = 0.62
        val a = 0.17
        val b = -0.13

        val result = OkColor().oklabToLinearSrgb(L, a, b)
        val expected = Triple(0.5228926220518537, 0.07239547193989224, 0.6266086535228347) // These values are placeholders

        assertEquals(expected.first, result.first, 0.0001, "Red channel does not match")
        assertEquals(expected.second, result.second, 0.0001, "Green channel does not match")
        assertEquals(expected.third, result.third, 0.0001, "Blue channel does not match")
    }

    @Test
    fun testFindCusp() {
        val a = 0.6
        val b = 0.2

        val result = OkColor().findCusp(a, b)
        val expectedL = 0.8  // This value is a placeholder
        val expectedC = 0.4  // This value is a placeholder

        assertEquals(result,expectedL)
//        assertTrue(abs(result.second - expectedC) < 0.01, "C_cusp does not match")
    }

    @Test
    fun testFindGamutIntersection() {
        val a = 0.5
        val b = -0.5
        val L1 = 0.9
        val C1 = 0.8
        val L0 = 0.5
        val cusp = Pair(0.8, 0.6)  // Example cusp values

        val result = OkColor().findGamutIntersection(a, b, L1, C1, L0, cusp)
        val expected = 0.7499999999999998

        assertEquals(expected, result, 0.01, "Intersection t does not match")
    }

    @Test
    fun testGetSTMax() {
        val a = 0.3
        val b = -0.3
        val cusp = Pair(0.9, 0.7) // Example cusp values

        val result = OkColor().getSTMax(a, b, cusp)
        val expectedS = 0.7777777778 // Placeholder for C/L
        val expectedT = 7.0          // Placeholder for C/(1-L)

        assertEquals(expectedS, result.first, 0.0001, "S does not match")
        assertEquals(expectedT, result.second, 0.0001, "T does not match")
    }

    @Test
    fun testGetSTMid() {
        val a = 0.2
        val b = -0.1

        val result = OkColor().getSTMid(a, b)
        val expectedS = 0.26945375922508963
        val expectedT = 0.6901280149386917

        assertEquals(expectedS, result.first, 0.0001, "S does not match")
        assertEquals(expectedT, result.second, 0.0001, "T does not match")
    }

    @Test
    fun srgbTransferFunction() {

        assertEquals(0.012920000000000001,OkColor().srgbTransferFunction(0.001))
        assertEquals(0.06100854008843736,OkColor().srgbTransferFunction(0.005))
    }

    @Test
    fun srgbTransferFunction_inv() {
        assertEquals(0.012920000000000001,OkColor().srgbTransferFunction_inv(0.001))
        assertEquals(0.06100854008843736,OkColor().srgbTransferFunction_inv(0.005))
    }

    @Test
    fun testOkhslToSrgb() {
        val h = 0.0 // Hue
        val s = 0.5 // Saturation
        val l = 0.5 // Lightness

        val result = OkColor().okhslToSrgb(h, s, l)
        val expected = Triple(128, 128, 128) // Expected sRGB values

        assertEquals(170,result.first)
        assertEquals(89,result.second)
        assertEquals(116,result.third)

    }

    @Test
    fun testGetCs() {
        val L = 0.5
        val a = 0.2
        val b = -0.1

        val results = OkColor().getCs(L, a, b)
        // Expected values should be computed or validated against a reliable source.
        val expectedC0 = 0.8   // Example expected value
        val expectedCmid = 0.9  // Example expected value
        val expectedCmax = 1.0  // Example expected value

        assertTrue { abs(results[0] - expectedC0) < 0.01 }
        assertTrue { abs(results[1] - expectedCmid) < 0.01 }
        assertTrue { abs(results[2] - expectedCmax) < 0.01 }
    }

}