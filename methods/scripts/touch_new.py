import os
from openpyxl import Workbook
from lxml import etree
import re

def writeExcel(path, data, name):
    wb = Workbook()
    ws = wb.active
    ws.append(['文件名', '序列名称', 'Region', '预测的功能类型', '序列开头', '序列结尾', '最相似', '相似的功能类型', '相似度', '氨基酸序列', '核酸序列'])
    for da in data:
        temp = [name]
        temp.extend(da)
        ws.append(temp)
    wb.save(path)

def handleH(path):
    try:
        with open(path, 'r') as p:
            page = etree.HTML(p.read())
        titles = page.xpath('//div[@class="record-overview-header"]/strong/text()')
        tbody = page.xpath('//div[@class="record-overview-details"]//table[@class="region-table"]//tbody')
        result = {}
        for title, tb in zip(titles, tbody):
            trs = tb.xpath('.//tr')
            tempp = []
            for tr in trs:
                tds = tr.xpath(".//td")
                temp = [title]
                for td in tds:
                    t = td.xpath(".//text()")
                    strs = ''.join(t).replace('&nbsp', " ").strip()
                    temp.append(strs)
                while len(temp) < 8:
                    temp.append('')
                tempp.append(temp)
            result[title] = tempp
        return result
    except Exception as e:
        print("HTML文件处理出错", e)

def handleF(path):
    try:
        with open(path, 'r') as p:
            page = p.read()
        cdss = re.findall(r'/translation="(.+?)"', page, flags=re.DOTALL)
        cds = ''.join([cd.replace("\n", "").replace(" ", "") for cd in cdss])

        origins = re.findall(r'ORIGIN(.+?)//', page, flags=re.DOTALL)[0]
        origin = re.findall(r'\d+? (.+?)\n', origins)
        origin = ''.join([cd.replace(" ", "") for cd in origin])
        return cds, origin
    except Exception as e:
        print("gbk文件处理出错", e)

def main():
    # 根目录文件夹名称
    root = "."

    fs = [i for i in os.listdir(root) if os.path.isdir(i)]
    for name in fs:
        folder = os.path.join(root, name)
        html = os.path.join(root, name, 'index.html')
        result = handleH(html)
        if not result:
            continue
        resultss = []
        for originalName in result:
            files = os.listdir(folder)
            for f in files:
                if originalName in f and "region" in f:
                    temp = []
                    index = int(re.findall(r'region(\d\d\d)', f)[0]) - 1
                    if 0 <= index < len(result[originalName]):
                        print("File:", f, "Index:", index)  # Debugging print statement
                        ffff = os.path.join(folder, f)
                        cds, origin = handleF(ffff)
                        temp.extend(result[originalName][index])
                        temp.append(cds)
                        temp.append(origin)
                        resultss.append(temp)
                    else:
                        print("Index out of range for file:", f)
        writeExcel("{}.xlsx".format(name), resultss, name)
        print(name, "finished")
    print("Done")

if __name__ == '__main__':
    main()
