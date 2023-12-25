const express = require('express')

const app = express();
const bodyParser = require('body-parser');
const cors = require('cors');
const multer = require('multer'); // npm package for handling file uploads
const {
    spawn
} = require('node:child_process');
const fs = require('fs')
const {
    parse
} = require("csv-parse");

app.use(bodyParser.urlencoded({
    extended: false
}));
app.use(bodyParser.json());
app.use(cors());

app.post('/building', () => {

})

app.get('/api', (req, res) => {
    return res.send('Hello World');
})

app.get('/api/get-response', async (req, res) => {
    let typeOfModel = fs.readFileSync('typeOfModel.txt', {
        encoding: 'utf8',
        flag: 'r'
    });
    let XZwater2 = [];
    fs.createReadStream("XZwater2.csv")
        .pipe(parse({
            delimiter: ",",
            from_line: 2
        }))
        .on("data", function (row) {
            row = [+row[0], +row[1] * -1]
            XZwater2.push(row);
        })

    let XZwater = [];
    fs.createReadStream("XZwater1.csv")
        .pipe(parse({
            delimiter: ",",
            from_line: 2
        }))
        .on("data", function (row) {
            row = [+row[0], +row[1] * -1]
            XZwater.push(row);
        })

    let nuWaterRight = [];
    fs.createReadStream("nuWaterRight.csv")
        .pipe(parse({
            delimiter: ",",
            from_line: 1
        }))
        .on("data", function (row) {
            row = row.map(k => +k);
            nuWaterRight.push(row);
        })
    let nuWaterLeft = [];
    fs.createReadStream("nuWaterLeft.csv")
        .pipe(parse({
            delimiter: ",",
            from_line: 1
        }))
        .on("data", function (row) {
            row = row.map(k => +k);
            nuWaterLeft.push(row);
        })
    let nuxy = [];
    fs.createReadStream("nuxy1.csv")
        .pipe(parse({
            delimiter: ",",
            from_line: 1
        }))
        .on("data", function (row) {
            row = row.map(k => +k);
            nuxy.push(row);
        })
    let nuxk = [];
    fs.createReadStream("nuxk.csv")
        .pipe(parse({
            delimiter: ",",
            from_line: 1
        }))
        .on("data", function (row) {
            row = row.map(k => +k);
            nuxk.push(row);
        })
    let Roka = [];
    fs.createReadStream("Roka1.csv")
        .pipe(parse({
            delimiter: ",",
            from_line: 2
        }))
        .on("data", function (row) {
            row = [+(+row[1]), +(+row[2])]
            Roka.push(row);
        })
    let XZsurface = [];
    fs.createReadStream("XZsurface.csv")
        .pipe(parse({
            delimiter: ",",
            from_line: 2
        }))
        .on("data", function (row) {
            row = [+(+row[1]), +(+row[2] * -1)]
            XZsurface.push(row);
        })
    let electrodes = [];
    return fs.createReadStream("Electrodes.csv")
        .pipe(parse({
            delimiter: ",",
            from_line: 2
        }))
        .on("data", function (row) {
            // console.log(row);
            row = [+(+row[1]), +(+row[2] * -1)]
            electrodes.push(row);
        }).on('end', () => {

            let response = {
                XZsurface,
                electrodes,
                Roka,
                nuxk,
                nuxy,
                nuWaterLeft,
                XZwater,
                nuWaterRight,
                XZwater2,
                typeOfModel,
            }
            return res.json(response)
        })

    // console.log(XZsurface)
})


app.post('/api/calculate', (req, res) => {
    const AbMetres = +req.body.AbMetres;
    const reliefStartPoint = +req.body.reliefStartPoint;
    const typeOfModel = +req.body.typeOfModel;
    const nks = +req.body.nks;
    const heightOfWaterAtLeftSide = +req.body.heightOfWaterAtLeftSide;
    const heightOfWaterAtRightSide = +req.body.heightOfWaterAtRightSide;
    let fortranCodes = spawn('/home/hadmin/diplomka-backend-master/a.out', []);
    let dataResult = [];
    let inputTimeout;

    const sendInput = (input) => {
        fortranCodes.stdin.write(input + '\n');
    };

    const setInputTimeout = (timeout, data) => {
        inputTimeout = setTimeout(() => {
            sendInput(data);
        }, timeout);
    };
    if ([1, 2].includes(typeOfModel)) {
        setInputTimeout(50, AbMetres);
        setInputTimeout(50, reliefStartPoint);
        setInputTimeout(50, typeOfModel);
        setInputTimeout(50, nks);
        setInputTimeout(50, heightOfWaterAtLeftSide);
    }
    if (typeOfModel === 3) {
        setInputTimeout(50, AbMetres);
        setInputTimeout(50, reliefStartPoint);
        setInputTimeout(50, typeOfModel);
        setInputTimeout(50, nks);
        setInputTimeout(50, heightOfWaterAtLeftSide);
        setInputTimeout(50, heightOfWaterAtRightSide);
    }

    fortranCodes.stdout.on('data', (data) => {
        console.log(`stdout: ${data}`);
        dataResult.push(data.toString());
    });

    fortranCodes.stderr.on('data', (data) => {
        console.error(`stderr: ${data}`);
    });

    console.log("IT WORKs")
    return fortranCodes.on('close', (code) => {
        console.log("IT WORKs ", code)

        // Handle exit code of the Fortran program
        console.log(`child process exited with code ${code}`);
        fs.writeFileSync('typeOfModel.txt', typeOfModel.toString());
        let dataForResponse = {
            data: dataResult,
            _success: true,
        }
        console.log("IT WORKs ui")

        return res.json(dataForResponse);
    });

})

app.get('/', (req, res) => {
    return res.send('Is worked')
})

const storage = multer.diskStorage({
    destination: function (req, file, cb) {
        cb(null, 'uploads/'); // specify the destination folder for uploaded files
    },
    filename: function (req, file, cb) {
        cb(null, file.fieldname + '-' + Date.now() + path.extname(file.originalname)); // generate unique filenames for uploaded files
    }
});
const upload = multer({
    storage
});

// Define a route for the file uploader endpoint
app.post('/upload', upload.single('file'), function (req, res, next) {
    if (!req.file) {
        return res.status(400).json({
            error: 'No file uploaded'
        });
    }

    // Access the uploaded file details from req.file object
    const {
        originalname,
        filename,
        path,
        size
    } = req.file;

    // Handle the uploaded file as needed (e.g., save to database, process, etc.)
    // ...

    // Send a response back to the client
    res.status(200).json({
        success: 'File uploaded successfully',
        file: req.file
    });
});



app.listen(3000, () => {
    console.log('hello world')
})