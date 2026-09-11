// QR helper for the poster, using only macOS's built-in Core Image - no downloads.
//
//   swift qr.swift gen "<text>"      prints the module matrix, one row per line ("1" = dark)
//   swift qr.swift read <image.png>  prints every QR payload Core Image finds in the image
import CoreImage
import Foundation

func generate(_ text: String) -> [[Bool]] {
    let filter = CIFilter(name: "CIQRCodeGenerator")!
    filter.setValue(text.data(using: .utf8), forKey: "inputMessage")
    filter.setValue("Q", forKey: "inputCorrectionLevel")   // ~25% damage tolerance, good for print
    let image = filter.outputImage!                          // one pixel per module
    let ctx = CIContext()
    let cg = ctx.createCGImage(image, from: image.extent)!
    let w = cg.width, h = cg.height
    var buf = [UInt8](repeating: 0, count: w * h)
    let gray = CGColorSpaceCreateDeviceGray()
    let bmp = CGContext(data: &buf, width: w, height: h, bitsPerComponent: 8, bytesPerRow: w,
                        space: gray, bitmapInfo: CGImageAlphaInfo.none.rawValue)!
    bmp.draw(cg, in: CGRect(x: 0, y: 0, width: w, height: h))
    var rows: [[Bool]] = []
    for y in 0..<h { rows.append((0..<w).map { buf[y * w + $0] < 128 }) }
    // Core Image pads the symbol with a light border; trim it so the caller controls the quiet zone.
    while let r = rows.first, !r.contains(true) { rows.removeFirst() }
    while let r = rows.last, !r.contains(true) { rows.removeLast() }
    let cols = rows[0].indices.filter { c in rows.contains { $0[c] } }
    return rows.map { Array($0[cols.first!...cols.last!]) }
}

func read(_ path: String) -> [String] {
    guard let img = CIImage(contentsOf: URL(fileURLWithPath: path)) else { return [] }
    let det = CIDetector(ofType: CIDetectorTypeQRCode, context: nil,
                         options: [CIDetectorAccuracy: CIDetectorAccuracyHigh])!
    return det.features(in: img).compactMap { ($0 as? CIQRCodeFeature)?.messageString }
}

let args = CommandLine.arguments
if args.count == 3 && args[1] == "gen" {
    for row in generate(args[2]) { print(row.map { $0 ? "1" : "0" }.joined()) }
} else if args.count == 3 && args[1] == "read" {
    for m in read(args[2]) { print(m) }
} else {
    FileHandle.standardError.write("usage: swift qr.swift gen <text> | read <image>\n".data(using: .utf8)!)
    exit(2)
}
