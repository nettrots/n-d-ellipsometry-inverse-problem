using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FieldsConnector
{
    public class FieldsList
    {

        private Dictionary<string, AbstractFieldConnector> _list;
        public AbstractFieldConnector this[string index]
        {

            get
            {
                return _list[index];

            }

        }
        public string indexOfInterface(object o)
        {
            foreach (var list in _list)
            {
                if (list.Value.checkForInterface(o)) return list.Key;
            }
            throw new Exception("No such interface bonded in fieldlist");
        }
        public void addBond(AbstractFieldConnector afc)
        {
            if(_list==null)_list=new Dictionary<string, AbstractFieldConnector>();
            _list.Add(afc.Name, afc);
        }
        public void SyncToAll()
        {
            foreach (var list in _list)
            {
                _list[list.Key].syncTo();
            }
        }
        public void SyncFromAll()
        {
            foreach (var list in _list)
            {
                _list[list.Key].syncFrom();
            }
        }

    }
}
